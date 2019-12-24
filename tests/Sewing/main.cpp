#include <cstring>
#include <iostream>
#include "nets.h"
#include "precomputation.h"
#include "computation.h"
#include "save-data.h"
#include <stdlib.h>
#include <math.h>
#include <stdio.h>
#include <sstream>
#include <vector>
#include <set>
#include <algorithm>
#include <gsl/gsl_multimin.h>
#include "Eigen/Core"
#include <autodiff/reverse.hpp>
#include <autodiff/reverse/eigen.hpp>
#include <chrono>
#include "include/MinEnergyDeformator.h"

extern "C" {
#include "libfrtprm.h"
#include "libaft.h"
#include "det.h"
}

#define ALPHA 1.0
#define BETA 0.7
#define SCALE 10

class Timer
{
private:
    using clock_t = std::chrono::high_resolution_clock;
    using second_t = std::chrono::duration<double, std::ratio<1> >;

    std::chrono::time_point<clock_t> m_beg;

public:
    Timer() : m_beg(clock_t::now())
    {
    }

    void reset()
    {
        m_beg = clock_t::now();
    }

    double elapsed() const
    {
        return std::chrono::duration_cast<second_t>(clock_t::now() - m_beg).count();
    }
};

double sew_fixed_mesh(net_t mesh);

using namespace autodiff;
using namespace Eigen;

class MeshInserter{
public:
    MeshInserter(net_t& mesh): mesh{mesh} {
        set_N();
    }

    struct MinimizationParams{
        double step_sz = 1.e-4;
        double tol = 1.e-4;
        double epsabs = 1.e-2;
        int maxits = 500;
    };

    struct SewEnergyParams{
        double alpha = 1.0;
        double beta = 1.0;
        double scale = 10.0;
        bool convexity = false;
    };


    net_t& mesh;
    std::vector<int> varMap;
    int varN = 0;
    SewEnergyParams energyParams;
    MinimizationParams minimizationParams;

    int insert(int freq = 1){
        size_t iter = 0;
        int status;

        const gsl_multimin_fdfminimizer_type *T;
        gsl_multimin_fdfminimizer *s;

        gsl_vector *x = gsl_vector_alloc (varN);
        for (int j = 0, i = 0; j < mesh.vrtx.count; ++j){
            if (is_fix(mesh.vrtx.nodes[j]->state)){
                continue;
            }
            for (int k = 0; k < 3; ++k)
                gsl_vector_set(x, i++, mesh.vrtx.nodes[j]->coord.coord[k]);// + (k == 2) ? 1.0 : 0.0);
        }

        gsl_multimin_function_fdf energy_func;

        energy_func.n = varN;
        energy_func.f = _f4gsl;
        energy_func.df = _df4gsl;
        energy_func.fdf = _fdf4gsl;
        energy_func.params = this;

        T = gsl_multimin_fdfminimizer_conjugate_fr;
        s = gsl_multimin_fdfminimizer_alloc (T, varN);

        gsl_multimin_fdfminimizer_set (s, &energy_func, x, minimizationParams.step_sz, minimizationParams.tol);

        do
        {
            iter++;
            status = gsl_multimin_fdfminimizer_iterate (s);

            if (status)
                break;

            status = gsl_multimin_test_gradient (s->gradient, minimizationParams.epsabs);

            if (status == GSL_SUCCESS)
                printf ("Minimum found at:\n");

            if ((freq > 0) && !(iter % freq))
                printf ("%5d: f() = %10.5f\n", iter, s->f);

        }
        while (status == GSL_CONTINUE && iter < minimizationParams.maxits);

        gsl_multimin_fdfminimizer_free (s);
        gsl_vector_free (x);

        return status;
    }

private:
    void set_N(){
        int i = 0;
        varMap.reserve(mesh.vrtx.count);
        for (int j = 0; j < mesh.vrtx.count; ++j){
            if (is_fix(mesh.vrtx.nodes[j]->state)){
                varMap[j] = -1;
                continue;
            }
            varMap[j] = i;
            i++;
        }
        varN = 3 * i;
    }

    static void set_coords_to_mesh(const gsl_vector* X, MeshInserter* m){
        net_t& mesh = m->mesh;
        for (int j = 0, i = 0; j < mesh.vrtx.count; ++j){
            if (is_fix(mesh.vrtx.nodes[j]->state)){
                continue;
            }
            for (int k = 0; k < 3; ++k)
                mesh.vrtx.nodes[j]->coord.coord[k] = gsl_vector_get(X, i++);
        }
    }



    static void _fdf4gsl(const gsl_vector *x, void *params, double *f, gsl_vector *df)
    {
        MeshInserter* m = (MeshInserter*)params;
        set_coords_to_mesh(x, m);
        *f = compute_energy(df, m);
    }

    static double _f4gsl (const gsl_vector *x, void *params)
    {
        MeshInserter* m = (MeshInserter*)params;
        set_coords_to_mesh(x, m);
        return compute_energy(NULL, m);
    }

    static void _df4gsl (const gsl_vector *x, void *params, gsl_vector *df)
    {
        MeshInserter* m = (MeshInserter*)params;
        set_coords_to_mesh(x, m);
        compute_energy(df, m);
    }

    static double compute_energy(gsl_vector * df, MeshInserter* m){
        if (df) gsl_vector_set_all(df, 0);

        net_t& mesh = m->mesh;
        std::vector<int>& varMap = m->varMap;
        double energy = 0;
        for (int i = 0; i < mesh.springs.count; ++i) {
            spring_t &sp = *mesh.springs.springs[i];
            if (is_fix(sp.ends[0]->state) && is_fix(sp.ends[1]->state)) continue;
            point_t p[2] = {sp.ends[0]->coord, sp.ends[1]->coord};
            double edge_grad[6];
            energy += loc_edge_energy(p, sp.l_0, edge_grad, m->energyParams, df != NULL);
            int emap[2] = {varMap[sp.ends[0]->id], varMap[sp.ends[1]->id]};
            if (sp.isdigedral) {
                point_t p[4];
                p[0] = sp.dihedral[0]->coord, p[1] = sp.dihedral[1]->coord;
                p[2] = sp.ends[0]->coord, p[3] = sp.ends[1]->coord;
                int dimap[4] = {varMap[sp.dihedral[0]->id], varMap[sp.dihedral[1]->id], varMap[sp.ends[0]->id],
                                varMap[sp.ends[1]->id]};
                double digedral_grad[12];
                energy += loc_digedral_energy(p, digedral_grad, m->energyParams, df != NULL);
                if (df != nullptr) {
                    for (int j = 0; j < 6; ++j)
                        edge_grad[j] += digedral_grad[j + 6];
                    for (int j = 0; j < 2; ++j) {
                        if (dimap[j] >= 0)
                            for (int k = 0; k < 3; ++k)
                                gsl_vector_set(df, dimap[j] * 3 + k,
                                               gsl_vector_get(df, dimap[j] * 3 + k) +
                                               digedral_grad[j * 3 + k]);
                    }
                }
            }
            if (df != nullptr) {
                for (int j = 0; j < 2; ++j) {
                    if (emap[j] >= 0)
                        for (int k = 0; k < 3; ++k)
                            gsl_vector_set(df, emap[j] * 3 + k,
                                           gsl_vector_get(df, emap[j] * 3 + k) +
                                           edge_grad[j * 3 + k]);
            }
            }
        }

        return energy;
    }

    static var coerced_digedral_energy(const VectorXvar& x, SewEnergyParams& ep)
    {
#define SS(i) x[0 + (0 + i)%3]*x[3 + (1 + i)%3]*x[6 + (2 + i)%3]
#define MM(i) x[0 + (5 - i)%3]*x[3 + (4 - i)%3]*x[6 + (3 - i)%3]
        auto abc = SS(0) + SS(1) + SS(2) - MM(0) - MM(1) - MM(2);
#undef SS
#undef MM
        auto a_2 = x[0] * x[0] + x[1] * x[1] + x[2] * x[2];
        auto a = sqrt(a_2);
#define DD(i, j) x[3*i + 0]*x[3*j + 0] + x[3*i + 1]*x[3*j + 1] + x[3*i + 2]*x[3*j + 2]
        auto bc = DD(1, 2);
        auto ac = DD(0, 2);
        auto ab = DD(0, 1);
#undef DD

        auto phi = atan(abc / ((bc - ab * ac / a_2) * a)) / M_PI;
        double beta = ep.beta;
        double scale = ep.scale;
        bool sign = ep.convexity;
        double k;
        if (sign)
            k = 0.1 * ((atan2(abc.get()->val, ((bc - ab * ac / a_2) * a).get()->val) > 0)? 1.0 : scale);
        else
            k = 0.1 * ((atan2(abc.get()->val, ((bc - ab * ac / a_2) * a).get()->val) < 0)? 1.0 : scale);

        return k * beta * phi * phi;
    }

    static double loc_digedral_energy(point_t p[4], double x[12], SewEnergyParams& ep, bool needDerivative){ //x - input-output variable
        int mask[3] = {3, 0, 1};

        VectorXvar xx(9);
        for (int k = 0; k < 3; ++k)
            for (int i = 0; i < 3; ++i)
                xx[3 * k + i] = p[mask[k]].coord[i] - p[2].coord[i];

        var y = coerced_digedral_energy(xx, ep);
        if (needDerivative) {
            VectorXd grad = gradient(y, xx);

            for (int i = 0; i < 3; ++i) {
                x[3 * 0 + i] = grad[3 * 1 + i];
                x[3 * 1 + i] = grad[3 * 2 + i];
                x[3 * 2 + i] = -(grad[3 * 0 + i] + grad[3 * 1 + i] + grad[3 * 2 + i]);
                x[3 * 3 + i] = grad[3 * 0 + i];
            }
        }

        return val(y);
    }

    static inline var __edge_energy(VectorXvar& y, SewEnergyParams& ep, double l0){
#define _ADD(I) (y[I] - y[I + 3]) * (y[I] - y[I + 3])
        var l = sqrt(_ADD(0) + _ADD(1) + _ADD(2));
#undef _ADD
        var err = (l - l0) / l0;
        double alpha = ep.alpha;
        var energy = alpha * err * err;
        return energy;
    }

    static double loc_edge_energy(point_t p[2], double l0, double x[6], SewEnergyParams& ep, bool needDerivative){
        VectorXvar y(6);
        for (int i = 0; i < 6; ++i)
            y[i] = p[i/3].coord[i%3];
        var energy = __edge_energy(y, ep, l0);
        if (needDerivative) {
            VectorXd grad = gradient(energy, y);
            for (int i = 0; i < 6; ++i) {
                x[i] = grad[i];
            }
        }

        return val(energy);
    }


};





net_t g_mesh;
int g_N = 0;
std::vector<int> g_map;
double g_timer;
Timer g_tt;
double _H_MESH = 1.0;

net_t generate_ozaki_template(int template_num, double mesh_h);
void sew_fixed_boundary(net_t& mesh, double R);
void set_coords_to_mesh(const gsl_vector* X);
void set_N();
double compute();
double compute1();
void compute_df(gsl_vector *df);
void print_length_quality();

double func(const gsl_vector* X){
    set_coords_to_mesh(X);
    return compute();
}

double func1(const gsl_vector* X){
    set_coords_to_mesh(X);
    return compute1();
}

double f4gsl1 (const gsl_vector *v, void *params)
{
    //std::cout << "f" << std::endl;
    double f = func1(v);

    return f;
}

double f4gsl (const gsl_vector *v, void *params)
{
    //std::cout << "f" << std::endl;
    double f = func(v);

    return f;
}

void df4gsl (const gsl_vector *v, void *params, gsl_vector *df)
{
    //std::cout << "df" << std::endl;
    //set_coords_to_mesh(v);
    compute_df(df);
}

void fdf4gsl(const gsl_vector *x, void *params, double *f, gsl_vector *df)
{
    //std::cout << "fdf"<< std::endl;
    *f = f4gsl(x, params);
    df4gsl(x, params, df);
}

int find_minimum_energy_df(double step_sz = 1.e-4, double tol = 1.e-4, double epsabs = 1.e-2, int maxits = 50000, double time/*secs*/ = 150){
    size_t iter = 0;
    int status;

    const gsl_multimin_fdfminimizer_type *T;
    gsl_multimin_fdfminimizer *s;

    double* par = NULL;

    gsl_vector *x = gsl_vector_alloc (g_N);
    for (int j = 0, i = 0; j < g_mesh.vrtx.count; ++j){
        if (is_fix(g_mesh.vrtx.nodes[j]->state)){
            continue;
        }
        for (int k = 0; k < 3; ++k)
            gsl_vector_set(x, i++, g_mesh.vrtx.nodes[j]->coord.coord[k]);// + (k == 2) ? 1.0 : 0.0);
    }

    gsl_multimin_function_fdf my_func;

    my_func.n = g_N;
    my_func.f = f4gsl;
    my_func.df = df4gsl;
    my_func.fdf = fdf4gsl;
    my_func.params = par;

    T = gsl_multimin_fdfminimizer_conjugate_fr;
    s = gsl_multimin_fdfminimizer_alloc (T, g_N);

    gsl_multimin_fdfminimizer_set (s, &my_func, x, step_sz, tol);

    Timer t;
   do
    {
        iter++;
        status = gsl_multimin_fdfminimizer_iterate (s);

        if (status)
            break;

        status = gsl_multimin_test_gradient (s->gradient, epsabs);

        if (status == GSL_SUCCESS)
            printf ("Minimum found at:\n");

        if (!(iter % 1))
        printf ("%5d: f() = %10.5f\n", iter, s->f);

    }
    while (status == GSL_CONTINUE && iter < maxits && t.elapsed() < time);

    gsl_multimin_fdfminimizer_free (s);
    gsl_vector_free (x);

    return status;
}

int find_minimum_energy(double init_sz = 1000, double drop_sz = 1.e-2, int maxits = 100000){
    double* par = NULL;

    const gsl_multimin_fminimizer_type *T = gsl_multimin_fminimizer_nmsimplex2;//gsl_multimin_fminimizer_nmsimplex2rand;
    gsl_multimin_fminimizer *s = NULL;
    gsl_vector *ss, *x;
    gsl_multimin_function minex_func;

    size_t iter = 0;
    int status;
    double size;

    /* Starting point */
    x = gsl_vector_alloc (g_N);
    for (int j = 0, i = 0; j < g_mesh.vrtx.count; ++j){
        if (is_fix(g_mesh.vrtx.nodes[j]->state)){
            continue;
        }
        for (int k = 0; k < 3; ++k)
            gsl_vector_set(x, i++, g_mesh.vrtx.nodes[j]->coord.coord[k]);// + (k == 2) ? 1.0 : 0.0);
    }


    /* Set initial step sizes to 0.1 */
    ss = gsl_vector_alloc (g_N);
    gsl_vector_set_all (ss, init_sz);

    /* Initialize method and iterate */
    minex_func.n = g_N;
    minex_func.f = f4gsl1;
    minex_func.params = par;

    s = gsl_multimin_fminimizer_alloc (T, g_N);
    gsl_multimin_fminimizer_set (s, &minex_func, x, ss);

    do
    {
        iter++;
        status = gsl_multimin_fminimizer_iterate(s);

        if (status)
            break;

        size = gsl_multimin_fminimizer_size (s);
        status = gsl_multimin_test_size (size, drop_sz);

        if (status == GSL_SUCCESS)
        {
            printf ("converged to minimum at\n");
        }
        if (! (iter % 1000))
        printf ("%5d: f() = %7.3f size = %.3f\n", iter, s->fval, size);
    }
    while (status == GSL_CONTINUE && iter < maxits);

    gsl_vector_free(x);
    gsl_vector_free(ss);
    gsl_multimin_fminimizer_free (s);

    return status;
}

using namespace autodiff;
using namespace Eigen;
var test_f(const VectorXvar& x){
    return sqrt(x[0] * x[0]+ x[1] * x[1] + x[2] * x[2]);
}

void test_compute_df(int i){
    VectorXvar x(3);
    x[0] = i;
    x[1] = 1;
    x[2] = 1;


    var y = test_f(x);
    VectorXd grad = gradient(y, x);
    //std::cout << "y = " << y << std::endl;
    if (i == 1)
    std::cout << "grad = \n" << grad << std::endl;
}

int main(int argc, const char* argv[]) {

    struct InputParams{
        double mesh_h = 1.0;
        int ModelNum = 19;
        double R_surf = 8.0;
        static InputParams processMainArgs(int argc, const char* argv[]){
            InputParams ip = InputParams();
            for (int i = 1; i < argc; i++) {
                //Print help message and exit
                if (strcmp(argv[i], "-h") == 0 || strcmp(argv[i], "--help") == 0) {
                    helpMessage:
                    {
                        std::cout << "Help message: " << std::endl;
                        std::cout << "Command line options: " << std::endl;
                        std::cout << "-sz, --mesh_size <Step of generating eigth sphere part mesh>" << std::endl;
                        std::cout << "-m, --model <Model number>" << std::endl;
                        std::cout << "-R, --surface_radius <Radius of cylinder for covering it by fixed boundary>" << std::endl;
                        std::cout << "-h, --help <Write help message>" << std::endl;
                    }
                    exit(0);
                }
                if (strcmp(argv[i], "-sz") == 0 || strcmp(argv[i], "--mesh_size") == 0){
                    std::cout << "mesh size = " << argv[i + 1] << std::endl;
                    ip.mesh_h = atof(argv[++i]);
                    continue;
                }
                if (strcmp(argv[i], "-m") == 0 || strcmp(argv[i], "--model") == 0){
                    std::cout << "model number = " << argv[i + 1] << std::endl;
                    int model_num = atoi(argv[++i]);
                    std::set<int> possible_sizes = {13, 15, 17, 19, 21, 23, 25, 27, 29, 31, 33, 35};
                    if (!possible_sizes.count(model_num)){
                        std::cout << "Wrong model num\n";
                        std::cout << "Only the following models are available: ";
                        for(auto i: possible_sizes)
                            std::cout << i << " ";
                        std::cout << std::endl;
                        exit(0);
                    }
                    ip.ModelNum = model_num;
                    continue;
                }
                if (strcmp(argv[i], "-R") == 0 || strcmp(argv[i], "--surface_radius") == 0){
                    std::cout << "surface_radius = " << argv[i + 1] << std::endl;
                    ip.R_surf = atof(argv[++i]);
                    continue;
                }
            }
            return ip;
        }
    };

    //test_compute_df();
    //return 0;
    /*VectorXvar y(12);
    y[0] = 0, y[1] = 0, y[2] = 1;
    y[3] = 0, y[4] = 1, y[5] = 0;
    y[6] = 0, y[7] = 0, y[8] = 0;
    y[9] = 1, y[10] = 0, y[11] = 0;
    var z = ddigedral_energy(y);
    std::cout << z << "\n" << gradient(z, y);
    return 0;*/

    InputParams params = InputParams::processMainArgs(argc, argv);
    net_t mesh = generate_ozaki_template(params.ModelNum, params.mesh_h);
    _H_MESH = params.mesh_h;
    g_mesh = mesh;
    set_N();
    net_t_set_elems_neighbours(mesh);
    nets_t _nets = nets_t_get_net(1); _nets.nets[0] = mesh;
    nets_t_set_relax_state(_nets, point_t_get_point(1, 0, 0));
    to_stl(_nets, "../bin/preresult");
    std::cout <<  "Energy_plate = " << compute() << std::endl;

    sew_fixed_boundary(mesh, params.R_surf);
    to_stl(_nets, "../bin/fixed");
    std::cout <<  "Energy_initial = " << compute() << std::endl;
    /*std::vector<double> deriv;
    deriv.resize(3 * g_N);
    compute_df(deriv.data());
    for (auto& i: deriv)
        std::cout << i << " ";
    std::cout << std::endl;
    return 0;*/
    Timer t;
//   find_minimum_energy_df();
    MinEnergyDeformator deformator(mesh);
    set_default_length_constr(deformator, ALPHA);
    set_default_digedral_angle_constr(deformator, BETA, SCALE);
    point_t normal = NORM(VEC(-0.4887521089927858, 0.16566805036803003, 0.8565485818343056));
    point_t n2 = NORM(VEC(0.7196727241521752, 0.055280565525260526, 0.6921092610177922));
    double b2 = DOT(n2, VEC(4.749989204287048, -10.906683337965926, 7.553886966850468));
    //SCAL_S(-1, &normal);
    double b = DOT(normal, VEC(-3.8634846210479736, -17.77840232849121, 7.0771379470825195));
    set_plane_constr(deformator, 1.0, normal, b);
    set_plane_constr(deformator, 1.0, n2, b2);
    std::cout <<  "Energy_initial = " << deformator.getEnergy() << std::endl;
    //exit(0);
    deformator.find_minimum_energy_df(1.e-4, 1.e-2, 5.e-2 );
    std::cout << "Autodif time = " << g_timer << std::endl;
    //MeshInserter m{mesh};
    //m.insert();
    std::cout << "t = " << t.elapsed() << std::endl;
    t.reset();

    //find_minimum_energy(1.0);
    //for (int i = 0; i < 5; ++i) {
        //find_minimum_energy(100, 1e-3, 500000);
        std::stringstream ss{};
        ss << "";//i;
        to_stl(_nets, ("../bin/result_" + ss.str()).c_str());
        std::cout << "Energy_result = " << compute1() << std::endl;
    print_length_quality();
    //}

    return 0;
}

void print_length_quality(){
    int q[20] = {};
    for (int i = 0; i < g_mesh.springs.count; ++i) {
        spring_t &sp = *g_mesh.springs.springs[i];
        if (is_fix(sp.ends[0]->state) && is_fix(sp.ends[1]->state)) continue;
        double l = LEN(DIF(sp.ends[0]->coord, sp.ends[1]->coord));
        double l_0 = sp.l_0;
        double dif = fabs(l - l_0) / l_0;
        for (int k = 1; k <= 20; ++k)
            if (dif > 0.05 * k)
                q[k - 1]++;
    }

    std::cout << "Quality 5%: ";
    for (int i = 1; i <= 20; ++i)
        std::cout << q[i - 1] << " ";
    std::cout << std::endl;
}

void set_coords_to_mesh(const gsl_vector* X){
    for (int j = 0, i = 0; j < g_mesh.vrtx.count; ++j){
        if (is_fix(g_mesh.vrtx.nodes[j]->state)){
            continue;
        }
        for (int k = 0; k < 3; ++k)
            g_mesh.vrtx.nodes[j]->coord.coord[k] = gsl_vector_get(X, i++);
    }
}

void set_coords_to_mesh(double* X){
    for (int j = 0, i = 0; j < g_mesh.vrtx.count; ++j){
        if (is_fix(g_mesh.vrtx.nodes[j]->state)){
            continue;
        }
        for (int k = 0; k < 3; ++k)
            g_mesh.vrtx.nodes[j]->coord.coord[k] = X[i++];
    }
}

double compute1(){
    double f = 0;
    double alpha = ALPHA;
    double beta = BETA;
    for (int i = 0; i < g_mesh.springs.count; ++i){
        spring_t& sp = *g_mesh.springs.springs[i];
        if (is_fix(sp.ends[0]->state) && is_fix(sp.ends[1]->state)) continue;
        double l = LEN(DIF(sp.ends[0]->coord, sp.ends[1]->coord));
        double l_0 = sp.l_0;
        f += alpha * (l - l_0) * (l - l_0) / l_0 / l_0 * 0; //TODO: zero
        if (sp.isdigedral) {
            point_t p[4];
            p[0] = sp.dihedral[0]->coord, p[1] = sp.dihedral[1]->coord;
            p[2] = sp.ends[0]->coord, p[3] = sp.ends[1]->coord;
            point_t n[2];
            for (int j = 0; j < 2; ++j) {
                point_t pp = point_t_or_area(p[2], p[3], p[j]);
                double len = LEN(pp);
                if (len > _H_MESH * 1.e-6)
                    n[j] = SCAL(1.0/len, pp);
                else
                    n[j] = ZERO();
            }
            SCAL_S(-1, &n[1]);
            point_t res = CROSS(n[0], n[1]);
            double cos_phi = DOT(n[0], n[1]);
            double sin_phi = LEN(res);
            if (DOT(res, DIF(p[3], p[2])) >= 0)
                sin_phi *= -1;
            double phi = atan2(sin_phi, cos_phi) / M_PI;
            double k = 0.1 * ((phi > 0)? 1.0 : SCALE);
            f += beta * k * phi * phi;
        }
    }

    return f;
    //return compute();
}

double compute(){
    double f = 0;
    double alpha = ALPHA;
    double beta = BETA;
    for (int i = 0; i < g_mesh.springs.count; ++i){
        spring_t& sp = *g_mesh.springs.springs[i];
        if (is_fix(sp.ends[0]->state) && is_fix(sp.ends[1]->state)) continue;
        double l = LEN(DIF(sp.ends[0]->coord, sp.ends[1]->coord));
        double l_0 = sp.l_0;
        double eps = (l - l_0) / l_0;
        f += alpha * eps * eps;
        if (sp.isdigedral) {
            point_t p[4];
            p[0] = sp.dihedral[0]->coord, p[1] = sp.dihedral[1]->coord;
            p[2] = sp.ends[0]->coord, p[3] = sp.ends[1]->coord;
            point_t n[2];
            for (int j = 0; j < 2; ++j) {
                point_t pp = point_t_or_area(p[2], p[3], p[j]);
                double len = LEN(pp);
                if (len > _H_MESH * 1.e-6)
                    n[j] = SCAL(1.0/len, pp);
                else
                    n[j] = ZERO();
            }
            SCAL_S(-1, &n[1]);
            point_t res = CROSS(n[0], n[1]);
            double cos_phi = DOT(n[0], n[1]);
            double sin_phi = LEN(res);
            if (DOT(res, DIF(p[3], p[2])) >= 0)
                sin_phi *= -1;
            double phi = atan2(sin_phi, cos_phi) / M_PI;
            double k = 0.1 * ((phi > 0)? 1.0 : SCALE);
            f += beta * k * phi * phi;
        }
    }

    return f;
}

// The scalar function for which the gradient is needed
var f(const VectorXvar& x, double *f)
{
#define SS(i) x[0 + (0 + i)%3]*x[3 + (1 + i)%3]*x[6 + (2 + i)%3]
#define MM(i) x[0 + (5 - i)%3]*x[3 + (4 - i)%3]*x[6 + (3 - i)%3]
    auto abc = SS(0) + SS(1) + SS(2) - MM(0) - MM(1) - MM(2);
#undef SS
#undef MM
    auto a_2 = x[0] * x[0] + x[1] * x[1] + x[2] * x[2];
    auto a = sqrt(a_2);
#define DD(i, j) x[3*i + 0]*x[3*j + 0] + x[3*i + 1]*x[3*j + 1] + x[3*i + 2]*x[3*j + 2]
    auto bc = DD(1, 2);
    auto ac = DD(0, 2);
    auto ab = DD(0, 1);
#undef DD
    auto phi = atan(abc / ((bc - ab * ac / a_2) * a)) / M_PI;
    double phival = atan2(abc.get()->val, -((bc - ab * ac / a_2) * a).get()->val);
    if (f) *f = phival / M_PI;

    double k = 0.1 * (((phival) > 0)? 1.0 : SCALE);

    return k * phi;
}

double deriv_f1(const point_t p[3], point_t res[3], double sc[2]){
    //f = (a, b, c) * |a| / ([a, b], [a, c])
    //p[0] <-> a
    //p[1] <-> b
    //p[2] <-> c
    double a2 = SQR_LEN(p[0]);
    double a = sqrt(a2);
    double ab = DOT(p[0], p[1]), ac = DOT(p[0], p[2]), bc = DOT(p[1], p[2]);
    point_t axb = CROSS(p[0], p[1]), cxa = CROSS(p[2], p[0]), bxc = CROSS(p[1], p[2]);
    double abc = DOT(axb, p[2]), axb_axc = -DOT(axb, cxa);
    //std::array<point_t&, 3> abc_p = {bxc, cxa, axb};
    double phi = abc * a;
    //point_t abc_adna = SCAL(abc/a, p[0]);
    //for (int i = 1; i < 3; ++i)
    //    SCAL_S(a, &abc_p);
    //point_t phi_da = SCAL_SUM(a, bxc, abc/a, p[0]);
    std::array<point_t, 3> phi_p = {SCAL_SUM(a, bxc, abc/a, p[0]), SCAL(a, cxa), SCAL(a, axb)};
    double psi = axb_axc;
    point_t psi_da;
    for (int i = 0; i < 3; ++i)
        psi_da.coord[i] = 2 * bc * p[0].coord[i] - ac * p[1].coord[i] - ab * p[2].coord[i];
    point_t psi_p[3] = {psi_da, SCAL_SUM(a2, p[2], -ac, p[0]), SCAL_SUM(a2, p[1], -ab, p[0])};
    double f = phi / psi;
    sc[0] = phi, sc[1] = -psi;
    double d1_psi = 1.0/psi, phid2psi = -phi / psi / psi;
    for (int i = 0; i < 3; ++i){
        res[i] = SCAL_SUM(d1_psi, phi_p[i], phid2psi, psi_p[i]);
    }
    return f;
}

double deriv_f2(const point_t p[3], point_t res[3], double beta, double scal){
    double sc[2];
    double x = deriv_f1(p, res, sc);
    double f = atan2(sc[0], sc[1]) / M_PI;
    double coef_arctg = 1.0 / (1.0 + x * x);
    double k = (f > 0) ? 0.1 : (scal * 0.1);
    double coef = coef_arctg * k * beta / M_PI;
    for (int i = 0; i < 3; ++i)
        SCAL_S(coef, &res[i]);
    return f;
}

double deriv_f3(const point_t p[4], double x[12], double beta, double scal){
    int mask[3] = {3, 0, 1};
    point_t pp[3];
    for (int i = 0; i < 3; ++i)
        pp[i] = ADD(p[mask[i]], -1, p[2]);
    point_t grad[3];
    double f = deriv_f2(pp, grad, beta, scal);

    for (int i = 0; i < 3; ++i) {
        x[3 * 0 + i] = grad[1].coord[i];
        x[3 * 1 + i] = grad[2].coord[i];
        x[3 * 2 + i] = -(grad[0].coord[i] + grad[1].coord[i] + grad[2].coord[i]);
        x[3 * 3 + i] = grad[0].coord[i];;
    }

    return f;
}

double digedral_df(point_t p[4], double x[12]){ //x - input-output variable
    int mask[3] = {3, 0, 1};

    VectorXvar xx(9);
    for (int k = 0; k < 3; ++k)
        for (int i = 0; i < 3; ++i)
            xx[3 * k + i] = p[mask[k]].coord[i] - p[2].coord[i];

    /*static int _id = 0; //TODO: delete this
    if (fabs((double)xx[3 * 0 + 2]) + fabs((double)xx[3 + 2]) + fabs((double)xx[3 * 2 + 2]) > DBL_EPSILON)
        _id++;*/


    double phi = 0;
    var y = f(xx, &phi);
    VectorXd grad = gradient(y, xx);

    VectorXd res(12);
    for (int i = 0; i < 3; ++i) {
        res[3 * 0 + i] = grad[3 * 1 + i];
        res[3 * 1 + i] = grad[3 * 2 + i];
        res[3 * 2 + i] = -(grad[3 * 0 + i] + grad[3 * 1 + i] + grad[3 * 2 + i]);
        res[3 * 3 + i] = grad[3 * 0 + i];
    }

    for (int i = 0; i < 12; ++i)
        x[i] = res[i];

    return phi;
}

void compute_df(gsl_vector *df){
    //std::fill(df, df + 3 * g_N, 0);
    gsl_vector_set_all(df, 0);

    double alpha = ALPHA;
    double beta = BETA;
    for (int i = 0; i < g_mesh.springs.count; ++i) {
        spring_t& sp = *g_mesh.springs.springs[i];
        double l = LEN(DIF(sp.ends[0]->coord, sp.ends[1]->coord));
        double l_0 = sp.l_0;
        int id[4] = {g_map[sp.ends[0]->id], g_map[sp.ends[1]->id], -1, -1};
        double coef = 2 * alpha * (l - l_0) / l / l_0 / l_0;
        for (int ii = 0; ii < 2; ++ii)
            for (int k = 0; k < 3; ++k) {
                if (id[ii] >= 0)
                    gsl_vector_set(df, id[ii] * 3 + k,
                                   gsl_vector_get(df, id[ii] * 3 + k) +
                                   coef * (sp.ends[0]->coord.coord[k] - sp.ends[1]->coord.coord[k])
                                                       * ((ii == 0) ? 1.0 : -1.0));
                /*df[id[ii] * 3 + k] += coef * (sp.ends[0]->coord.coord[k] - sp.ends[1]->coord.coord[k])
                        * ((ii == 0) ? 1.0 : -1.0);*/
            }

        if (sp.isdigedral) {
            point_t p[4];
            p[0] = sp.dihedral[0]->coord; id[0] = g_map[sp.dihedral[0]->id];
            p[1] = sp.dihedral[1]->coord; id[1] = g_map[sp.dihedral[1]->id];
            p[2] = sp.ends[0]->coord;     id[2] = g_map[sp.ends[0]->id];
            p[3] = sp.ends[1]->coord;     id[3] = g_map[sp.ends[1]->id];
            double x[12];
            g_tt.reset();
            double phi = deriv_f3(p, x, 1.0, SCALE);//
            double y[12];
            /*double phi0 =  digedral_df(p, y);
            for (int i = 0; i < 12; ++i)        //TODO: comparator
                if (fabs(x[i] -y[i]) > sqrt(DBL_EPSILON)) {
                    std::cout << i << ": ERROR\n";
                    break;
                }*/
            g_timer += g_tt.elapsed();

            /*point_t n[2];
            for (int j = 0; j < 2; ++j) {
                point_t pp = point_t_or_area(p[2], p[3], p[j]);
                double len = LEN(pp);
                if (len > _H_MESH * 1.e-6)
                    n[j] = SCAL(1.0/len, pp);
                else
                    n[j] = ZERO();
            }
            SCAL_S(-1, &n[1]);
            point_t res = CROSS(n[0], n[1]);
            double cos_phi = DOT(n[0], n[1]);
            double sin_phi = LEN(res);
            if (DOT(res, DIF(p[3], p[2])) >= 0)
                sin_phi *= -1;
            double phi = atan2(sin_phi, cos_phi) / M_PI;
            if (fabs(phi - phi0) > sqrt(DBL_EPSILON))
                std::cout << "ERROR\n";*/
            double k = 1;//0.1 * ((phi > 0)? 1.0 : SCALE);
            double _edge_coef = 2 * beta * k * phi;

            for (int ii = 0; ii < 4; ++ii)
                if (id[ii] >= 0)
                for (int kk = 0; kk < 3; ++kk)
                    gsl_vector_set(df, id[ii] * 3 + kk,
                                   gsl_vector_get(df, id[ii] * 3 + kk) -
                                           _edge_coef * x[ii*3 + kk]);
                    //df[id[ii] * 3 + k] += 2 * beta * k * phi * x[ii*3 + k];
        }
    }

   /* for (int i = 0; i < g_N; ++i)
        std::cout << i << ": " << gsl_vector_get(df, i) << std::endl;*/
   //TODO: почему в направлении x всегда 0?

}

void set_N(){
    int i = 0;
    g_map.reserve(g_mesh.vrtx.count);
    for (int j = 0; j < g_mesh.vrtx.count; ++j){
        if (is_fix(g_mesh.vrtx.nodes[j]->state)){
            g_map[j] = -1;
            continue;
        }
        g_map[j] = i;
        i++;
    }
    g_N = 3 * i;
}



void sew_fixed_boundary(net_t& mesh, double R){
    double x_max = 0;
    for (int i = 0; i < mesh.vrtx.count; ++i){
        double x = mesh.vrtx.nodes[i]->coord.coord[0];
        if (x_max < fabs(x)) x_max = fabs(x);
    }
    if (x_max >= M_PI * R) {
        std::cout << "Too small cylinder radius. You should use R > " << x_max / M_PI << " for this model" << std::endl;
        exit(0);
    }
    for (int i = 0; i < mesh.vrtx.count; ++i){
        node_t& node = *mesh.vrtx.nodes[i];
        if (is_fix(node.state)){
            double phi = node.coord.coord[0] / R;
            node.coord.coord[0] = R * sin(phi);
            //node.coord.coord[1] /= 2;
            node.coord.coord[2] = R * cos(phi);
        }
        else {
            //node.coord.coord[1] += 15;
            //node.coord.coord[2] += 5;
            node.coord.coord[2] += 1.1*R;      //make artificial shift
        }
    }
}

//ozaki genereator functions
static double global_size = 0.1;

static struct {
    double a, b, alpha, beta, radius, shift;
} ozaki;

/* this function controls the desired size of the mesh in space */
/* x, y, z  -- coords of point in space                         */
static double fsize(double x, double y, double z) {
    return global_size;
    (void) x,  (void) y,  (void) z;
}

/* Surface parametric function */
static int surface_param(int i, double u, double v, double *px, double *py, double *pz) {
    *px = u,  *py = v,  *pz = 0.0;
    return 1;
    (void) i;
}
/* Line parametric function */
static int line_param(int i, double t, double *pu, double *pv) {
    double ax = 0, ay = ozaki.a * sin(ozaki.beta);
    double bx = ozaki.a * cos(ozaki.beta), by = 0;
    double cx = ozaki.a * cos(ozaki.beta) - ozaki.b*sin(ozaki.alpha), cy = -ozaki.b*cos(ozaki.alpha);
    double dx = -ozaki.a * cos(ozaki.beta) + ozaki.b*sin(ozaki.alpha), dy = -ozaki.b*cos(ozaki.alpha);
    double ex = -ozaki.a * cos(ozaki.beta), ey = 0;
    if (i == 1) {
        *pu = (1-t)*ax + t*bx,  *pv = (1-t)*ay + t*by;
    } else if (i == 2) {
        *pu = (1-t)*bx + t*cx,  *pv = (1-t)*by + t*cy;
    } else if (i == 3) {
        *pu = (1-t)*cx + t*dx,  *pv = (1-t)*cy + t*dy;
    } else if (i == 4) {
        *pu = (1-t)*dx + t*ex,  *pv = (1-t)*dy + t*ey;
    } else if (i == 5) {
        *pu = (1-t)*ex + t*ax,  *pv = (1-t)*ey + t*ay;
    } else if (i == 6) {
        *pu = ozaki.radius*cos(ozaki.alpha + t*(M_PI-2*ozaki.alpha)),  *pv = -ozaki.radius*sin(ozaki.alpha + t*(M_PI-2*ozaki.alpha)) + ozaki.shift;
    }
    return 1;
}

int create_ozaki(int *pnV, double *vertex, int *pnF, int *face, int *facematerial, int *pnE, int *edge, int *edgematerial, int nnV, int nnF, int nnE, int ozaki_size, double size) {
    int nVVert = 5,  nLine = 5,  nSurface = 1;
    double VVert[5*3] = {-60, -60, -1,  60, 60, 1,  0,0,0, 0,0,0, 0,0,0};
    int LineD[5*3] = {1, 2, 1,  2, 3, 1,  3, 4, 1,  4, 5, 1,  5, 1, 1}; /*edges*/
    int LineP[5*2] = {1, 1,  1, 2,  1, 6,  1, 4,  1, 5}; /*param functions*/
    int exportCurves[5] = {1, 2, 2, 2, 1}; /*curve colors*/
    double LineT[5*2] = {0, 1,  0, 1,  0, 1, 0, 1,  0, 1}; /*parameter pairs*/
    int SurfL[5] = {5 /*edges*/, 1 /*param*/, 1, 0 /*backcolor*/, 0 /*direction*/}; /*main surface*/
    int SurfI[5*2] = {1, 1 /*direction*/,  2, 1 /*direction*/,  3, 1,  4, 1,  5, 1};
    double SurfT[4] = {-60.0,  60.0, -60.0, 60.0}; /*parametrization bbox*/

    if (ozaki_size == 13) {
        ozaki.alpha = 5.0 * M_PI / 180.0,  ozaki.beta = 11.4 * M_PI / 180.0,  ozaki.a = 10.6,  ozaki.b = 11.0;
    } else if (ozaki_size == 15) {
        ozaki.alpha = 5.0 * M_PI / 180.0,  ozaki.beta = 10.9 * M_PI / 180.0,  ozaki.a = 11.5,  ozaki.b = 11.2;
    } else if (ozaki_size == 101) { //my size
        ozaki.alpha = 3.1 * M_PI / 180.0,  ozaki.beta = 10.9 * M_PI / 180.0,  ozaki.a = 13.1,  ozaki.b = 11;
    } else if (ozaki_size == 17) {
        ozaki.alpha = 5.0 * M_PI / 180.0,  ozaki.beta =  9.0 * M_PI / 180.0,  ozaki.a = 12.8,  ozaki.b = 13.3;
    } else if (ozaki_size == 19) {
        ozaki.alpha = 5.0 * M_PI / 180.0,  ozaki.beta =  8.8 * M_PI / 180.0,  ozaki.a = 13.8,  ozaki.b = 13.1;
    } else if (ozaki_size == 21) {
        ozaki.alpha = 5.0 * M_PI / 180.0,  ozaki.beta =  7.8 * M_PI / 180.0,  ozaki.a = 14.7,  ozaki.b = 13.6;
    } else if (ozaki_size == 23) {
        ozaki.alpha = 5.0 * M_PI / 180.0,  ozaki.beta =  7.7 * M_PI / 180.0,  ozaki.a = 15.8,  ozaki.b = 13.5;
    } else if (ozaki_size == 25) {
        ozaki.alpha = 5.0 * M_PI / 180.0,  ozaki.beta = 10.0 * M_PI / 180.0,  ozaki.a = 16.9,  ozaki.b = 13.5;
    } else if (ozaki_size == 27) {
        ozaki.alpha = 5.0 * M_PI / 180.0,  ozaki.beta = 10.0 * M_PI / 180.0,  ozaki.a = 17.9,  ozaki.b = 13.6;
    } else if (ozaki_size == 29) {
        ozaki.alpha = 5.0 * M_PI / 180.0,  ozaki.beta =  9.3 * M_PI / 180.0,  ozaki.a = 18.8,  ozaki.b = 13.7;
    } else if (ozaki_size == 31) {
        ozaki.alpha = 5.0 * M_PI / 180.0,  ozaki.beta =  8.3 * M_PI / 180.0,  ozaki.a = 19.9,  ozaki.b = 14.0;
    } else if (ozaki_size == 33) {
        ozaki.alpha = 5.0 * M_PI / 180.0,  ozaki.beta =  8.2 * M_PI / 180.0,  ozaki.a = 20.8,  ozaki.b = 14.1;
    } else if (ozaki_size == 35) {
        ozaki.alpha = 5.0 * M_PI / 180.0,  ozaki.beta =  8.0 * M_PI / 180.0,  ozaki.a = 21.8,  ozaki.b = 14.1;
    } else {
        fprintf(stderr, "unknown ozaki size");
    }
    ozaki.radius = (ozaki.a * cos(ozaki.beta) - ozaki.b*sin(ozaki.alpha)) / cos(ozaki.alpha);
    ozaki.shift = -ozaki.b*cos(ozaki.alpha) + ozaki.radius * sin(ozaki.alpha);

    global_size = size;
    return ani3d_surface_edges_boundary_(&nVVert, VVert, &nLine, LineD, LineP, LineT, &nSurface, SurfL, SurfI, SurfT,
                                         NULL, surface_param, line_param, NULL, fsize, exportCurves,
                                         pnV, vertex, pnF, face, facematerial, pnE, edge, edgematerial,
                                         &nnV, &nnF, &nnE
    );
}

static int mark_vertices(int nV, int *vertexmaterial, int nE, int *edge, int *edgematerial, int defcolor) {
    int i;
    for (i=0; i<nV; i++)  vertexmaterial[i] = defcolor;
    for (i=0; i<nE; i++) {
        vertexmaterial[edge[2*i+0]-1] |= edgematerial[i];
        vertexmaterial[edge[2*i+1]-1] |= edgematerial[i];
    }

    vertexmaterial[1] = 2;
    vertexmaterial[4] = 2;

    return 0;
}


static double size = 1.0;

static int write_custom_format(char *fname, int nV, double *vertex, int *vertexmaterial, int nF, int *face) {
    FILE *f = fopen(fname, "w");
    int i;
    double n[3], q;
    if (!f)  return perror(fname),  -1;
    fprintf(f, "%d %d\n", nV, nF);
    for (i=0; i<nV; i++) {
        fprintf(f, "%le %le %le  %d\n", vertex[3*i+0], vertex[3*i+1], vertex[3*i+2], vertexmaterial[i]);
    }
    for (i=0; i<nF; i++) {
        normvec3i(vertex, face[3*i+0]-1, face[3*i+1]-1, face[3*i+2]-1, n);
        q = sqrt(n[0]*n[0] + n[1]*n[1] + n[2]*n[2]);
        if (q > 0)  n[0] /= q,  n[1] /= q,  n[2] /= q;
        fprintf(f, "%d %d %d  %le %le %le\n", face[3*i+0], face[3*i+1], face[3*i+2], n[0], n[1], n[2]);
    }
    fclose(f);
    return 0;
}

net_t generate_net_from_ani_mesh(int nV, double* vertex, int* vertexmaterial, int nF, int* face, int nE, int* edge){
    vrtx_t vrtx = vrtx_t_construct(nV);
    int mask[3] = {IN, MOB_BND, FIX_BND};
    for (int i = 0; i < nV; i++){
        int state = mask[vertexmaterial[i]];
        double h = 0.5;
        vrtx.nodes[i] = node_t_construct(vertex[3*i+0], vertex[3*i+1], vertex[3*i+2], h, state, i);
    }
    elems_t elems = elems_t_construct(nF);
    for (int i = 0; i < nF; i++){
        elems.elems[i] = elem_t_construct(vrtx.nodes[face[3*i+0]-1], vrtx.nodes[face[3*i+2]-1], vrtx.nodes[face[3*i+1]-1], i);
    }
    springs_t sprs = springs_t_construct((3 * nF - nE) / 2 + nE);
    net_t net = net_t_get(vrtx, elems, sprs);
    net_t_set_springs(&net);
    init_net(&net);

    return net;
}

/* generator*/
net_t generate_ozaki_template(int template_num, double mesh_h) {
    int    nnF = 600000, nnV = 600000, nnE = 600000;
    int     nF = 0,        nV = 0,        nE = 0;
    int    *face  = 0, *facematerial = 0;
    int    *edge  = 0, *edgematerial = 0;
    double *vertex = 0;
    int    *vertexmaterial = 0;

    int ozaki_size = template_num;
    size = mesh_h;

    // allocate memory for mesh structures
    vertex         = (double*)malloc(sizeof(double) * 3 * nnV);
    vertexmaterial = (int*)   malloc(sizeof(int)        * nnV);
    face           = (int*)   malloc(sizeof(int)    * 3 * nnF);
    facematerial   = (int*)   malloc(sizeof(int)        * nnF);
    edge           = (int*)   malloc(sizeof(int)    * 4 * nnE);
    edgematerial   = (int*)   malloc(sizeof(int)        * nnE);

    printf("size = %lf\n", size);

    printf("\n * Generating surface mesh\n"); /* just an example */
    create_ozaki(&nV, vertex, &nF, face, facematerial, &nE, edge, edgematerial,  nnV, nnF, nnE,  ozaki_size, size);
    mark_vertices(nV, vertexmaterial, nE, edge, edgematerial, 0);

    // Write surface mesh
    write_front_gmv   ("../bin/surface.gmv", nV, vertex, nF, face, facematerial);
    write_front       ("../bin/surface.smv", nV, vertex, nF, face, facematerial); // for smv
    write_custom_format("../bin/surface.txt", nV, vertex, vertexmaterial, nF, face);

    printf("INFO: nV = %d, nF = %d, nE = %d\n", nV, nF, nE);

    net_t net = generate_net_from_ani_mesh(nV, vertex, vertexmaterial, nF, face, nE, edge);
    free(vertex),  free(vertexmaterial),  free(face),  free(facematerial),  free(edge),  free(edgematerial);
    return net;
}

