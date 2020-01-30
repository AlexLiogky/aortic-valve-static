#include <cstring>
#include <iostream>
#include "nets.h"
#include "precomputation.h"
#include "computation.h"
#include "save-data.h"
#include <stdlib.h>
#include <math.h>
#include <stdio.h>
#include <vector>
#include <functional>
#include <string>
#include <fstream>
#include <gsl/gsl_poly.h>
#include <float.h>
#include "nlohmann/json.hpp"

using namespace std;

extern "C" {
#include "libfrtprm.h"
#include "libaft.h"
#include "det.h"
}

//#define  RENDERING
//#define _RADIUS_TESTER
#include <SDL/SDL.h>
#include <GL/gl.h>
#include <GL/glu.h>
#include "camera.h"


extern "C" {
unsigned int get_Cnt_Nodes();
point_t get_Mid_Shift(int i);
void restart_Mid_Shift();
void update_Mid_Shift(nets_t nets);
void construct_Mid_Shift(nets_t nets);
double get_Mid_Shift_len(int nsteps);
}

net_t generate_eigth_sphere_part(double R0, double h, std::vector<int>& vm, string prename);
net_t generate_cilinder_part(double R0, double l, double phi, double size, std::vector<int>& vm, string prename);
void set_thickness(net_t& mesh, double thickness);

void test_compute_nets(nets_t nets, double P, double delta, unsigned int n, double eps, int freq, std::vector<int>& vrt_mask, function<void(node_t*, int)> next_updater);
double test_spheriable(net_t mesh); //return approximate radius of sphere
double test_cilindrable(net_t mesh);//return approximate radius of cilinder
double evaluate_radius(int meshtype, net_t mesh);
double predicted_radius(double P /*Hg mm*/, double R0 /*mm*/, double T /*Pa*/, ElasticModels mod, int mesh_type,  double* prms);

struct Lendl{};
struct Logger{
    ostream& out = cout;
    ofstream log;
    bool log_b = false;

    void open( const std::string &filename,
               ios_base::openmode mode = ios_base::out ){
            log.open(filename, mode);
            log << scientific;
            log.precision(14);
            log_b = true;
    }
    void close(){
        if (log_b)
            log.close();
        log_b = false;
    }
    template <typename T>
    Logger& operator<< (const T& value){
        if (log_b)
            log << value;
        out << value;
        return *this;
    }
    template <typename T>
    void to_log (const T& value){
        if (log_b)
            log << value;
    }
    Logger& operator<< (const Lendl& value){
        if (log_b)
            log << endl;
        out << endl;
        return *this;
    }
    Logger& operator<< (const point_t& p){
        operator<<("point(");
        for (int i = 0; i < 2; ++i)
            *this << (p.coord[i]) << ", ";
        *this << (p.coord[2]) << ")";
    }
};
Lendl lendl;
Logger Log;

struct InputParams{
    using json = nlohmann::json;

    double R0 = 10;
    double l = 10;
    double phi = M_PI_2;
    double h_mesh = 0.4;
    int MeshType = 1;
    ElasticModels EModelType = EMOD_REINFORCING;
    double thickness = 0.5;
    double eprms[10] = {1000000, 2.3, 0, 0};
    double rel_allow_shift = 0.02 * 4;
    double rel_max_shift = 0.05 * 4;
    double P = 80; //mm Hg
    double delta = 5e-7;
    double err = 1e-5;
    double maxits = 100000;
    int compute_print_freq = 100;
    string dir = "../bin/";
    string prefix = "gent_1_";
    string logfile = "logfile.txt";
    json init_data;
private:
    static json get_json(string filename){
        std::ifstream js;
        js.open(filename);
        if (! js.good()) {
            throw std::runtime_error("Unsuccessful attempt to open file \"" + filename + "\"");
        }
        json conf;
        js >> conf;
        js.close();

        return std::move(conf);
    }
    static ElasticModels getModel(string ModelName){
        transform(ModelName.begin(), ModelName.end(), ModelName.begin(), ::tolower);
        map<string, ElasticModels> models = {
                {"msm", EMOD_MSM},
                {"mass-spring", EMOD_MSM},
                {"mass-spring model", EMOD_MSM},
                {"trqs", EMOD_TRQS},
                {"st venan", EMOD_TRQS},
                {"neogook", EMOD_NEOGOOK},
                {"neogokean", EMOD_NEOGOOK},
                {"neo-gookean", EMOD_NEOGOOK},
                {"reinforcing", EMOD_REINFORCING},
                {"gent", EMOD_REINFORCING},
                {"compressible neogook", EMOD_REINFORCING}
        };
        if (!models.count(ModelName)) return EMOD_MSM;
        else return models[ModelName];
    }
    void processConfigFile(const char* filename){
        json conf = get_json(filename);
        R0 = conf["R"];
        l = conf["l"];
        double _phi = conf["phi"];
        phi = _phi / 180 * M_PI;
        h_mesh = conf["mesh_step"];
        MeshType = conf["mesh_type"];
        EModelType = getModel(conf["elastic_model"]);
        thickness = conf["thickness"];
        vector<double> v = conf["model_params"].get<vector<double >>();
        for (int i = 0; i < 10 && i < v.size(); ++i)
            eprms[i] = v[i];
        rel_allow_shift = conf["rel_allow_shift"];
        rel_max_shift = conf["rel_max_shift"];
        P = conf["P"];
        delta = conf["delta"];
        err = conf["velocity_decay"];
        maxits = conf["max_count_of_iterations"];
        compute_print_freq = conf["frequency_of_printing"];
        dir = conf["save_directory"];
        prefix = conf["file_prefix"];
        logfile = conf["logfile_name"];
        init_data = std::move(conf);
    }
public:
    static InputParams processMainArgs(int argc, const char* argv[]){
        InputParams ip = InputParams();
        int i = 0;
        if (argc == 1) goto helpMessage;
        for (i = 1; i < argc; i++) {
            //Print help message and exit
            if (strcmp(argv[i], "-h") == 0 || strcmp(argv[i], "--help") == 0) {
                helpMessage:
                std::cout << "Help message: " << std::endl;
                std::cout << "Command line options: " << std::endl;
                std::cout << "-c, --config <Main configuration parameters file name>" << std::endl;
                exit(0);
            }
            if (strcmp(argv[i], "-c") == 0 || strcmp(argv[i], "--config") == 0) {
                std::cout << "Main configuration file found: " << argv[i + 1] << std::endl;
                ip.processConfigFile(argv[i + 1]);
                i++;
                continue;
            }
        }
        return ip;
    }
};

function<void(node_t*, int)> get_next_updater(InputParams& prms){
    if (prms.MeshType == 0){
        return [](node_t* n, int mask){
            for (int p = 0; p < 3; ++p) {
                if (!mask) break;
                if (mask % 2) n->next.coord[p] = 0;
                mask /= 2;
            }
        };
    } else if (prms.MeshType == 1){
        auto conv = [phi = prms.phi](node_t* n, int mask){
            if (mask % 2){
                n->next.coord[1] = n->coord.coord[1];//0
            }
            mask /= 2;
            if (mask % 2){
                n->next.coord[2] = n->coord.coord[2];
            }
            mask /= 2;
            if (mask % 2){
                point_t norm = VEC(-sin(phi), cos(phi), 0);
                point_t shift = DIF(n->next, n->initial);
                ADD_S(&(n->next), -DOT(shift, norm), norm);
            }
        };
        return conv;
    }else
        return [](node_t* n, int mask){};
}

int generate_config_files(string dir);

int main(int argc, const char* argv[]) {
    InputParams params = InputParams::processMainArgs(argc, argv);
    //return generate_config_files(params.dir);
    if (params.logfile != "") Log.open(params.dir + params.prefix + params.logfile);
//    double _max = 2.0 / 1000;
//    for (int i = 0; i <= 1000; ++i){
//        double kk = _max * i;
//        //double P = kk * (params.eprms[0] / 3) * params.thickness / (params.R0 * 133.3);
//        Log << kk << ";";
//    }
//    Log << lendl;
//    for (int i = 0; i <= 1000; ++i) {
//        double kk = _max * i;
//        double P = kk * (params.eprms[0] / 3) * params.thickness / (params.R0 * 133.3);
//        double R = predicted_radius(P, params.R0, params.thickness, params.EModelType, params.MeshType, params.eprms);
//        Log << R << ";";
//    }
//    exit(0);
    Log.to_log(params.init_data); Log << lendl;
    std::vector<int> vrtmaterial;
    net_t mesh;
    if (params.MeshType == 0)
        mesh = generate_eigth_sphere_part(params.R0, params.h_mesh, vrtmaterial, params.dir + params.prefix);
    else if (params.MeshType == 1)
        mesh = generate_cilinder_part(params.R0, params.l, params.phi, params.h_mesh, vrtmaterial, params.dir + params.prefix);
    else {
        Log << "ERROR: MeshType = " << params.MeshType << " is not defined" << lendl;
        return -1;
    }
    set_thickness(mesh, params.thickness);

    net_t_set_elems_neighbours(mesh);		//order is unimportant
    nets_t _nets = nets_t_get_net(1); _nets.nets[0] = mesh;
    set_relax_consts(_nets, params.rel_allow_shift, params.rel_max_shift);
    nets_t_set_relax_state(_nets, point_t_get_point(1, 0, 0));
    to_stl(_nets, (params.dir + params.prefix + "preresult").c_str());
    //exit(0);
    set_default_elastic_model(params.EModelType);
    set_default_elastic_params(params.eprms, 10);
    double R_exact = predicted_radius(params.P, params.R0, params.thickness, params.EModelType, params.MeshType, params.eprms);
    double kk = 133.3 * params.P * params.R0 / (params.eprms[0] / 3 * params.thickness);
    Log << "P * R / ( mu * H) = " << kk << lendl;
    double R_init = evaluate_radius(params.MeshType, mesh);
    test_compute_nets(_nets, params.P, params.delta, params.maxits, params.err, params.compute_print_freq, vrtmaterial, get_next_updater(params));
    to_stl(_nets, (params.dir + params.prefix + "result").c_str());
    double R_model = evaluate_radius(params.MeshType, mesh);
    Log << "Evaluate radius = " << R_model << lendl;
    Log << "Exact radius = " << R_exact << lendl;
    Log << "|R_exact - R_eval| = " << fabs(R_exact - R_model) << lendl;

    return 0;
}

double solve_poly_equation(const int n, double* a){
    gsl_poly_complex_workspace * w = gsl_poly_complex_workspace_alloc (n);
    double z[2 * (n - 1)];
    gsl_poly_complex_solve (a, n, w, z);
    gsl_poly_complex_workspace_free (w);
    int i = 0;
    int real = 0;
    double rreal = z[1];
    double eps = 1.0e-7;
    int solved = 0;
    for (i = 0; i < n - 1; ++i)
        if ((0 < z[2*i] && z[2*i] < 1 + eps)) {
            real = i;
            rreal = z[2*i + 1];
            if (fabs(rreal) < eps)
                solved = 1;
            break;
        }

    for (; i < n - 1; ++i)
        if ((0 < z[2*i] && z[2*i] < 1 + eps) && (fabs(z[2*i+1]) < fabs(rreal)) && (fabs(z[2*i+1]) < eps)){
            real = i;
            rreal = z[2*i + 1];
            solved = 2;
        }
    if (!solved || fabs(z[real * 2 + 1]) > eps)
#ifndef _RADIUS_TESTER
        throw std::runtime_error("Failed to solve the equation");
#else
        return 0;
#endif
    return (z[2*real] < 1) ? z[2*real] : 1;
}

double pred_lambda_neogook_cilinder(double pr_hmu){
    if (pr_hmu >= 1)
#ifndef _RADIUS_TESTER
        throw std::runtime_error("Failed to solve the equation, cilinder will grow to infininity");
#else
        return 0.0;
#endif
    return 1 / pow(1 - pr_hmu, 0.25);
}

double pred_lambda_gent_sphere(double pr_hmu, double Jm){
    double k = pr_hmu / 2;
    const int n = 10;
    double a[n] = {2*k, 0, -k*(Jm + 3), Jm, 0, 0, k, 0, 0, -Jm};
    return 1/solve_poly_equation(n, a);
}

double pred_lambda_gent_cilinder(double pr_hmu, double Jm){
    double k = pr_hmu;
    if (fabs(k) < DBL_EPSILON) return 1.0;
    double Jm_k = Jm / k;
    double a = Jm_k - (Jm + 2), b = 1, c = -Jm_k;
    array<double, 3> x;
    int n = gsl_poly_solve_cubic(a, b, c, &x[0], &x[1], &x[2]);
    sort(x.begin(), x.begin() + n);
    double z = 0; bool finded = false;
    for (int i = 0; i < n; ++i)
        if (x[i] > 1) {
            z = x[i];
            finded = true;
            break;
        }
    if (!finded)
#ifndef _RADIUS_TESTER
        throw std::runtime_error("Failed to solve the equation, there are no root > 1");
#else
        return 0.0;
#endif

    return sqrt(z);
}

double pred_lambda_neogook_sphere(double pr_hmu){
    double k = pr_hmu / 2;
    const int n = 8;
    double a[n] = { -k, 1, 0, 0, 0, 0, 0, -1 };

    return 1/solve_poly_equation(n, a);
}

double predicted_radius(double P /*Hg mm*/, double R0 /*mm*/, double T /*Pa*/, ElasticModels mod, int mesh_type,  double* prms){
    double mu = prms[0] / 3;
    P *= 133.3;
    double k = P * R0 / (mu * T);
    if (mod == EMOD_NEOGOOK && mesh_type == 0){
        return R0 * pred_lambda_neogook_sphere(k);
    }
    if (mod == EMOD_NEOGOOK && mesh_type == 1){
        return R0 * pred_lambda_neogook_cilinder(k);
    }
    if (mod == EMOD_REINFORCING && mesh_type == 0){
        return R0 * pred_lambda_gent_sphere(k, prms[1]);
    }
    if (mod == EMOD_REINFORCING && mesh_type == 1){
        return R0 * pred_lambda_gent_cilinder(k, prms[1]);
    }
    return -1e20;
}

double evaluate_radius(int meshtype, net_t mesh){
    if (meshtype == 0) return test_spheriable(mesh);
    if (meshtype == 1) return test_cilindrable(mesh);
    return 0;
}

double test_spheriable(net_t mesh){
    double r = 0;
    for (int i = 0; i < mesh.vrtx.count; ++i)
    {
        r += LEN(mesh.vrtx.nodes[i]->coord);
    }
    r /= mesh.vrtx.count;
    double mid_div = 0, max_div = 0;
    for (int i = 0; i < mesh.vrtx.count; ++i)
    {
        double div = fabs(LEN(mesh.vrtx.nodes[i]->coord) - r);
        mid_div += div;
        max_div = (max_div >= div) ? max_div : div;
    }
    mid_div /= mesh.vrtx.count;
    Log << "Evaluate radius = " << r << lendl;
    Log << "Maximal deviation from radius is " << max_div << lendl;
    Log << "Middle deviation from radius is " << mid_div << lendl;

    return r;
}

double test_cilindrable(net_t mesh){
    double r = 0;
    auto length = [](const point_t& p){
        return sqrt(p.coord[0] * p.coord[0] + p.coord[1] * p.coord[1]);
    };
    for (int i = 0; i < mesh.vrtx.count; ++i)
    {
        r += length(mesh.vrtx.nodes[i]->coord);
    }
    r /= mesh.vrtx.count;
    double mid_div = 0, max_div = 0;
    for (int i = 0; i < mesh.vrtx.count; ++i)
    {
        double div = fabs(length(mesh.vrtx.nodes[i]->coord) - r);
        mid_div += div;
        max_div = (max_div >= div) ? max_div : div;
    }
    mid_div /= mesh.vrtx.count;
    Log << "Evaluate radius = " << r << lendl;
    Log << "Maximal deviation from radius is " << max_div << lendl;
    Log << "Middle deviation from radius is " << mid_div << lendl;

    return r;
}

#ifdef RENDERING
point_t double_to_rgb(double x, double scale)
{
    x /= scale;
    point_t p = {};
    if (x < 0)
    {
        p.coord[2] = 1;
        return p;
    }
    if (x < 0.5)
    {
        p.coord[2] = 1 - 2* x;
        p.coord[1] = 2 * x;
        return p;
    }
    if (x < 1)
    {
        p.coord[1] = 1 - 2* (x - 0.5);
        p.coord[0] = 2 * (x - 0.5);
        return p;
    }
    p.coord[0] = 1;
    return p;
}

point_t get_color(const net_t& net, node_t* node)
{
    double r = 0;
    for (int i = 0, cnt = node->cnt_springs; i < cnt; ++i)
    {
        spring_t& spr = *net.springs.springs[node->springs_id[i]];
        r += fabs(spr.l - spr.l_0) / spr.l_0;
    }
    r /= node->cnt_springs;
    const double max_deform = 0.6;
    return double_to_rgb(r, max_deform);
}

void renderSoftBody(const net_t net, point_t fcolor){
    glColor3d(fcolor.coord[0], fcolor.coord[1], fcolor.coord[2]);
    glBegin(GL_TRIANGLES);
    for (unsigned int i = 0; i < net.elems.count; ++i)
    {
        elem_t* elem = net.elems.elems[i];
        for(int j = 0; j < 3; ++j)
        {
            fcolor = get_color(net, elem->vrts[j]);
            glColor3d(fcolor.coord[0], fcolor.coord[1], fcolor.coord[2]);
            point_t x = elem->vrts[j]->coord;
            glVertex3d(x.coord[0], x.coord[1], x.coord[2]);
        }
    }
    glEnd();

    glColor3f(0.6, 0.6, 0.6);
    glBegin(GL_LINES);
    for (unsigned int i = 0; i < net.springs.count; ++i)
    {
        spring_t* spr = net.springs.springs[i];
        for(int j = 0; j < 2; ++j)
        {
            point_t x = spr->ends[j]->coord;
            glVertex3d(x.coord[0], x.coord[1], x.coord[2]);
        }
    }
    glEnd();
}

void display(nets_t nets, camera_t* cam)
{
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
    glLoadIdentity();
    camera_t_Control(cam);
    //drawSkybox(50);
    camera_t_UpdateCamera(cam);
    renderSoftBody(nets.nets[0], point_t_get_point(0.296, 0.221, 0.231));
}
#endif

void test_compute_nets(nets_t nets, double P, double delta, unsigned int n, double eps, int freq, std::vector<int>& vrt_mask, function<void(node_t*, int)> next_updater){
    construct_Mid_Shift(nets);
    unsigned int i, k = 0;
    double init_div = 0;
    double res = 0;
    int crush = 0;


    bool running = true;
#ifdef RENDERING
    SDL_Init(SDL_INIT_EVERYTHING);
    SDL_SetVideoMode(1280,800,32,SDL_OPENGL);
    Uint32 _start;
    SDL_Event event;
    float angle=45;
    glClearColor(0,0,0,1);
    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
    gluPerspective(angle,1280.0/800.0,1,1000);
    glMatrixMode(GL_MODELVIEW);
    glEnable(GL_DEPTH_TEST);
    camera_t* cam = camera_t_construct(point_t_get_point(0.800262, -2.22123, 3.87111));//-45.8, -22, 34.6));
    cam->camPitch = 33.2;//-8;
    cam->camYaw = 9;//308;
#endif

    for (i = 1; i <= n && k < 3 && running; i++){
        int flag = i % 1000;
#ifdef RENDERING
        _start=SDL_GetTicks();
        while(SDL_PollEvent(&event))
        {
            switch(event.type)
            {
                case SDL_QUIT:
                    running=0;
                    break;
                case SDL_KEYDOWN:
                    switch(event.key.keysym.sym)
                    {
                        case SDLK_ESCAPE:
                            running=0;
                            break;
                        case SDLK_y:
                            camera_t_mouseIn(cam, 0);
                            break;
                        case SDLK_SPACE:
                        {
                            Log << "current possition: " << cam->loc << lendl;
                            Log << "camPitch = " << cam->camPitch << ", camYaw = " << cam->camYaw << lendl;
                            break;
                        }
                        default: break;
                    }
                    break;
                case SDL_MOUSEBUTTONDOWN:
                {
                    camera_t_mouseIn(cam, !camera_t_isMouseIn(cam));
                    break;
                }
                default: break;

            }
        }
#endif
        {
            {
                unsigned int nets_cnt = nets.count, ii;
                double shift = 0;
                int ret = 0;
                for (ii = 0; ii < nets_cnt; ii++) {
                    if (net_is_static(nets.nets[ii])) continue;
                    double max_net_shift = compute_free_nexts(nets.nets[ii], P, delta);
                    if (shift < max_net_shift) shift = max_net_shift;
                }

                if (shift > Max_shift && FLAG_REL) {
                    relaxation(nets, sqrt(Allow_shift / shift));
                    ret++;
                }
                crush += ret;
            }

            for (int j = 0; j < nets.nets[0].vrtx.count; ++j) {
                int mask = vrt_mask[j];
                next_updater(nets.nets[0].vrtx.nodes[j], mask);
            }

            if (!flag) {
                double diviation = get_diviation(nets);
                Log << "it = " << i << ", diviation = " << diviation / delta << ", relation = " << diviation / init_div << lendl;
            }

            update_Mid_Shift(nets);

            if (i == 1) {
                double max_Mid_Shift = 0;
                unsigned int l;
                for (l = 0; l < get_Cnt_Nodes(); l++) {
                    double len = SQR_LEN(get_Mid_Shift(l));
                    if (max_Mid_Shift < len) max_Mid_Shift = len;
                }
                Log << "Max_Shift = " << sqrt(max_Mid_Shift) << lendl;
                init_div = get_diviation(nets);
            }
            if (!(i % freq)) {
                res = get_Mid_Shift_len(freq);
                double rms = res / get_Cnt_Nodes(); // delta ;
                Log << "RMS shift per iter of node = " << rms << " mm, relation = " << res / init_div << ", cr = " << crush << lendl;
                crush = 0;
                restart_Mid_Shift();
                if (res < eps * init_div) k++;
                else k = 0;
            }

            update_nets(nets);
        }
#ifdef RENDERING
        display(nets, cam);
        SDL_GL_SwapBuffers();
        /*if(1000.0/60>SDL_GetTicks()-_start)
            SDL_Delay(1000.0/60-(SDL_GetTicks()-_start));*/
#endif
    }
    Log << "Made " << --i << " iterations, relation = " << res / init_div << lendl;

}

void set_thickness(net_t& mesh, double thickness){
    for (int i = 0; i < mesh.vrtx.count; ++i)
        mesh.vrtx.nodes[i]->h = thickness;
}


//mesh generator functions//
static double global_radius = 1.0;
static double global_size = 0.1;
static double global_len = 1.0;
static double global_phi = M_PI_2;

/* this function controls the desired size of the mesh in space */
/* x, y, z  -- coords of point in space                         */
static double fsize(double x, double y, double z) {
    return global_size;
    (void) x,  (void) y,  (void) z;
}

/* Surface parametric function */
static int surface_param(int i, double u, double v, double *px, double *py, double *pz) {
    double z;
    z = 1.0 - u*u - v*v;
    if (z < 0.0)  z = 0.0;
    z = global_radius * sqrt(z);
    *px = global_radius * u,  *py = global_radius * v,  *pz = ((i % 2) ? -1 : 1) * z;
    return 1;
}
static double periodic(int i, int d) {
    return  0.0;
    (void) i,  (void) d;
}
/* Line parametric function */
static int line_param(int i, double t, double *pu, double *pv) {
    if (i == 1) {
        *pu = 0.0,  *pv = sin(t);
    } else if (i == 3) {
        *pu = cos(t),  *pv = sin(t);
    } else if (i == 2){
        *pu = cos(t), *pv = 0;
    }
    return 1;
}

int create_sphere_quarter(int *pnV, double *vertex, int *pnF, int *face, int *facematerial, int *pnE, int *edge, int *edgematerial, int nnV, int nnF, int nnE, double radius, double size) {

    int nVVert = 3,  nLine = 3,  nSurface = 1;
    double VVert[] = {-2*radius, -2*radius, -2*radius,  2*radius, 2*radius, 2*radius,  0, 0, 2*radius};  /* these will be recomputed, just define the bbox */
    int LineD[] = {2, 3, 1,  3, 1, 1,  1, 2, 1}; /*edges*/
    int LineP[] = {1, 1,  1, 2,  1, 3}; /*param functions*/
    int exportCurves[] = {1, 2, 3}; /*curve colors*/
    double LineT[] = {M_PI/2, 0,  M_PI/2, 0,  0, M_PI/2 }; /*parameter pairs*/
    int SurfL[] = {3 /*edges*/, 1 /*param*/, 1, 0 /*backcolor*/, 0 /*direction*/}; /*main surface*/
    int SurfI[] = {1, 0 /*direction*/,  2, 0 /*direction*/,  3, 0 /*direction*/};
    double SurfT[4] = {-1.0,  1.0, -1.0, 1.0}; /*parametrization bbox*/

    global_radius = radius;
    global_size = size;
    return ani3d_surface_edges_boundary_(&nVVert, VVert, &nLine, LineD, LineP, LineT, &nSurface, SurfL, SurfI, SurfT,
                                         NULL, surface_param, line_param, periodic, fsize, exportCurves,
                                         pnV, vertex, pnF, face, facematerial, pnE, edge, edgematerial,
                                         &nnV, &nnF, &nnE
    );
}

static int mark_vertices(int nV, int *vertexmaterial, int nE, int *edge, int *edgematerial, int defcolor) {
    int i;
    for (i=0; i<nV; i++)
        vertexmaterial[i] = defcolor;
    for (i=0; i<nE; i++) {
        vertexmaterial[edge[2*i+0]-1] |= edgematerial[i];
        vertexmaterial[edge[2*i+1]-1] |= edgematerial[i];
    }
    return 0;
}

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
    for (int i = 0; i < nV; i++){
        int state = (vertexmaterial[i] == 0) ? IN : MOB_BND;
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

net_t generate_eigth_sphere_part(double R0, double size, std::vector<int>& vm, string prename){
    Log << "Generating of eighth sphere part" << lendl;
    int    nnF = 600000, nnV = 600000, nnE = 600000;
    int     nF = 0,        nV = 0,        nE = 0;
    int    *face  = 0, *facematerial = 0;
    int    *edge  = 0, *edgematerial = 0;
    double *vertex = 0;
    int    *vertexmaterial = 0;

    // allocate memory for mesh structures
    vertex         = (double*)malloc(sizeof(double) * 3 * nnV);
    vertexmaterial = (int*)   malloc(sizeof(int)        * nnV);
    face           = (int*)   malloc(sizeof(int)    * 3 * nnF);
    facematerial   = (int*)   malloc(sizeof(int)        * nnF);
    edge           = (int*)   malloc(sizeof(int)    * 4 * nnE);
    edgematerial   = (int*)   malloc(sizeof(int)        * nnE);

    Log << "\tMesh step = " << size << lendl;

    Log << "\n\t* Generating surface mesh" << lendl;
    create_sphere_quarter(&nV, vertex, &nF, face, facematerial, &nE, edge, edgematerial,  nnV, nnF, nnE,  R0, size);
    mark_vertices(nV, vertexmaterial, nE, edge, edgematerial, 0);
    // Write surface mesh
    write_front_gmv    ((char*)(prename + "eigth_sphere.gmv").c_str(), nV, vertex, nF, face, facematerial);
    write_front        ((char*)(prename + "eigth_sphere.smv").c_str(), nV, vertex, nF, face, facematerial); // for smv
    write_custom_format((char*)(prename + "eigth_sphere.txt").c_str(), nV, vertex, vertexmaterial, nF, face);

    Log << "\t\"INFO: nV = " << nV << ", nF = " << nF << ", nE = " << nE << lendl;

    /* Add layers */
    Log << "\n\t* Adding layers" << lendl;

    vm.resize(nV);
    for (int i = 0; i < 3; ++i){
        double err = size * 2e-7;
        int mask[4] = {0, 1, 2, 4};
        vm[i] = 0;
        for (int j = 0; j < 3; ++j)
            if (fabs(vertex[3*i + j]) < err) vm[i] |= mask[j+1];
    }
    for (int i = 3; i < nV; ++i){
        int mask[4] = {0, 1, 2, 4};
        vm[i] = mask[vertexmaterial[i]];
    }

    net_t net = generate_net_from_ani_mesh(nV, vertex, vertexmaterial, nF, face, nE, edge);
    free(vertex),  free(vertexmaterial),  free(face),  free(facematerial),  free(edge),  free(edgematerial);
    return net;
}

/* Surface parametric function */
static int surface_param_1(int i, double u, double v, double *px, double *py, double *pz) {
    *px = global_radius * cos(u);
    *py = global_radius * sin(u);
    *pz = v;
    //cout << "surf: " << *px << " " << *py << " " << *pz << " (" << u << ", "<< v << ")" << endl;

    return 1;
}

/* Line parametric function */
static int line_param_1(int i, double t, double *pu, double *pv) {
    if (i == 1) {
        *pu = 0.0,  *pv = t;
    } else if (i == 3) {
        *pu = global_phi,  *pv = t;
    } else if (i == 2){
        *pu = t, *pv = global_len/2;
    }
    else if (i == 4){
        *pu = t, *pv = -global_len/2;
    }
    //cout << "line: " << i << " " << *pu << " " << *pv << endl;
    return 1;
}

int create_cilinder_part(int *pnV, double *vertex, int *pnF, int *face, int *facematerial, int *pnE, int *edge, int *edgematerial, int nnV, int nnF, int nnE,
        double radius, double len, double phi, double size) {

    int nVVert = 4,  nLine = 4,  nSurface = 1;
    double bb = len + radius;
    double r = radius, l = len;
    double VVert[] = {r, 0, -l/2,  r, 0, l/2,  r*cos(phi), r*sin(phi), l/2,  r*cos(phi), r*sin(phi), -l/2};
            //{-2*bb, -2*bb, -2*bb,  2*bb, 2*bb, 2*bb,  0, 0, 0,  bb, bb, -bb};  /* these will be recomputed, just define the bbox */

    int LineD[] = {1, 2, 1,  2, 3, 1,  3, 4, 1,  4, 1, 1}; /*edges*/
    int LineP[] = {1, 1,  1, 2,  1, 3,  1, 4}; /*param functions*/
    int exportCurves[] = {1, 2, 4, 2}; /*curve colors*/
    double LineT[] = {-l/2, l/2,  0, phi,  l/2, -l/2,  phi, 0}; /*parameter pairs*/

    int SurfL[] = {4 /*edges*/, 1 /*param*/, 1, 0 /*backcolor*/, 1 /*direction*/}; /*main surface*/
    int SurfI[] = {1, 0 /*direction*/,  2, 0 /*direction*/,  3, 0 /*direction*/,  4, 0};
    double SurfT[4] = {0,  phi, -l/2, l/2}; /*parametrization bbox*/

    global_radius = radius;
    global_size = size;
    global_len = len;
    global_phi = phi;
    return ani3d_surface_edges_boundary_(&nVVert, VVert, &nLine, LineD, LineP, LineT, &nSurface, SurfL, SurfI, SurfT,
                                         NULL, surface_param_1, line_param_1, periodic, fsize, exportCurves,
                                         pnV, vertex, pnF, face, facematerial, pnE, edge, edgematerial,
                                         &nnV, &nnF, &nnE
    );
    //    LineD[3*nLine+0] = v1; LineD[3*nLine+1] = v2; // indices of start and end points
//    LineD[3*nLine+2] = 2; // this curve belongs to two (2) surfaces
// first one have parametrization function number 1 (F_1)
// and in this surface curve could be parametrized by V_5(u), where u are from -1.0 to 1.0
//LineP[2*m+0] = 1; LineP[2*m+1] = 5; LineT[2*m+0] = -1.0; LineT[2*m+1] =  1.0; m++;
// Now we will define surfaces
// For each surface we define 5 integer numbers (SurfL):
// 1st: number of boundary curves
// 2nd: number of parametrization function
// 3rd: color of the face
// 4th: 0 if the face is boundary face, or the second color of the face if it is used to split volume
// 5th: 0 if orientation of the face corresponds with parametrization, 1 if it should be inverted, 0 for flat faces
// For each surface we alse define minimax values of (u,v) parametrization (SurfT)
// u_min, u_max, v_min, v_max
// For each boundary curve we define a pair of integers (SurfI):
// 1st: index of curve
// 2nd: orientation,  0 for normal, 1 for reversed
//	Here is an example how we can add surface:
/*
	SurfL[5*nSurface+0] = 2; // surface is bounded by two curves
	SurfL[5*nSurface+1] = 1; // surface is parametrized by F_1
	SurfL[5*nSurface+2] = 1; // the color of face is 1
	SurfL[5*nSurface+3] = 0; // the face is boundary
	SurfL[5*nSurface+4] = 0; // orientation of the face corresponds with parametrization
	// minimax values of parametrs
	SurfT[4*nSurface+0] = -1.0; SurfT[4*nSurface+1] = -1.0; SurfT[4*nSurface+2] =  1.0; SurfT[4*nSurface+3] =  1.0;
	// first curve is curve number c1, orientation is normal
	SurfI[2*m+0] = c1; SurfI[2*m+1] = 0; m++;
	// second curve is curve number c2, orientation is reversed
	SurfI[2*m+0] = c2; SurfI[2*m+1] = 1; m++;
	nSurface++;
*/
// In common case to invert orientation of the face one should invert the 5th parameter in SurfL
// and invert orientation of all boundary curves
// To simplify our program we can use the following macros for flat quadrilaterals:
// Here for each curve (V,I) is for index of the curve V and orientation I
}

net_t generate_cilinder_part(double R0, double l, double phi, double size, std::vector<int>& vm, string prename){
    Log << "Generating of cilinder part with R = " << R0 << ", l = " << l << ", phi = " << phi << lendl;
    int    nnF = 600000, nnV = 600000, nnE = 600000;
    int     nF = 0,        nV = 0,        nE = 0;
    int    *face  = 0, *facematerial = 0;
    int    *edge  = 0, *edgematerial = 0;
    double *vertex = 0;
    int    *vertexmaterial = 0;

    // allocate memory for mesh structures
    vertex         = (double*)malloc(sizeof(double) * 3 * nnV);
    vertexmaterial = (int*)   malloc(sizeof(int)        * nnV);
    face           = (int*)   malloc(sizeof(int)    * 3 * nnF);
    facematerial   = (int*)   malloc(sizeof(int)        * nnF);
    edge           = (int*)   malloc(sizeof(int)    * 4 * nnE);
    edgematerial   = (int*)   malloc(sizeof(int)        * nnE);

    Log << "\tMesh step = " << size << lendl;

    Log << "\n\t* Generating surface mesh" << lendl;
    create_cilinder_part(&nV, vertex, &nF, face, facematerial, &nE, edge, edgematerial,  nnV, nnF, nnE,  R0, l, phi, size);
    mark_vertices(nV, vertexmaterial, nE, edge, edgematerial, 0);
    // Write surface mesh
    write_front_gmv    ((char*)(prename + "cilinder_part.gmv").c_str(), nV, vertex, nF, face, facematerial);
    write_front        ((char*)(prename + "cilinder_part.smv").c_str(), nV, vertex, nF, face, facematerial); // for smv
    write_custom_format((char*)(prename + "cilinder_part.txt").c_str(), nV, vertex, vertexmaterial, nF, face);

    Log << "\t\"INFO: nV = " << nV << ", nF = " << nF << ", nE = " << nE << lendl;

    /* Add layers */
    Log << "\n\t* Adding layers" << lendl;

    vm.resize(nV);
    copy(vertexmaterial, vertexmaterial + nV, vm.begin());

    net_t net = generate_net_from_ani_mesh(nV, vertex, vertexmaterial, nF, face, nE, edge);
    free(vertex),  free(vertexmaterial),  free(face),  free(facematerial),  free(edge),  free(edgematerial);
    return net;
}

int generate_config_files(string dir){
    map<pair<int, ElasticModels>, double> brds;
    brds.insert({{0, EMOD_NEOGOOK}, 0.6});
    brds.insert({{1, EMOD_NEOGOOK}, 0.99});
    brds.insert({{0, EMOD_REINFORCING}, 1.98});
    brds.insert({{1, EMOD_REINFORCING}, 1.98});
    array<double, 4> mesh_steps = {0.1, 0.2, 0.4, 0.8};
    array<int, 2> mesh_types = {0, 1};
    array<pair<ElasticModels, string>, 2> e_models = {pair{EMOD_NEOGOOK, "neogook"}, pair{EMOD_REINFORCING, "gent"}};
    vector<tuple<double, int, string, double>> variants;
    variants.reserve(4 * 2 * 2 * 10);
    for(auto& i: e_models)
        for (auto& j: mesh_types){
            double brd = brds[{j, i.first}];
            for (int k = 1; k <= 10; ++k)
                for (auto l: mesh_steps)
                    variants.push_back(tuple{l, j, i.second, (brd * k) / 10});
        }
    for (int i = 0; i < variants.size(); ++i){
        //if ((i + 1) % 4 != 1) continue;
        auto& v = variants[i];
        using json = nlohmann::json;
        json js;
        js["R"] = 10;
        js["l"] = 7;
        js["phi"] = 90;
        js["mesh_step"] = get<0>(v);
        js["mesh_type"] = get<1>(v);
        js["elastic_model"] = get<2>(v);
        js["thickness"] = 0.5;
        js["model_params"] = {1e6, 2.3};
        js["rel_allow_shift"] = 0.08;
        js["rel_max_shift"] = 0.2;
        double P = get<3>(v) * 1e6 / 3 * 0.5 / (10 * 133.3);
        js["P"] = P;
        double _P = ((int)(P * 100)) / 100;
        js["delta"] = 5e-7;
        js["velocity_decay"] = 1e-5;
        js["max_count_of_iterations"] = 40000000;
        js["frequency_of_printing"] = 250;
        js["save_directory"] = "../bin/";
        js["file_prefix"] = to_string(get<1>(v)) + "_" + get<2>(v) + "_" + to_string((int)_P*100) + "_" + to_string((int)(get<0>(v)*10));
        js["logfile_name"] = "logfile.txt";
        std::ofstream ff;
        ff.open(dir + "config_" + to_string(i+1) + ".json");
        ff << js;
        ff.close();
    }
    return 0;
}