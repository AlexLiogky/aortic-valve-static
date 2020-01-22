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
#include <gsl/gsl_poly.h>

extern "C" {
#include "libfrtprm.h"
#include "libaft.h"
#include "det.h"
}

//#define  RENDERING
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

net_t generate_eigth_sphere_part(double R0, double h, std::vector<int>& vm);
void set_thickness(net_t& mesh, double thickness);

void test_compute_nets(nets_t nets, double P, double delta, unsigned int n, double eps, int freq, std::vector<int>& vrt_mask);
double test_spheriable(net_t mesh); //return approximate radius
double predicted_radius(double P, double R0, double T);

int main(int argc, const char* argv[]) {

    struct InputParams{
        double R0 = 10;
        double h_mesh = 0.4;
        int ModelNum = -1;
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
                        std::cout << "-R0, -R0= <Radius of sphere part>" << std::endl;
                        std::cout << "-m, --model <Model number>" << std::endl;
                        std::cout << "-h, --help <Write help message>" << std::endl;
                    }
                    exit(0);
                }
                if (strcmp(argv[i], "-sz") == 0 || strcmp(argv[i], "--mesh_size") == 0){
                    std::cout << "mesh size = " << argv[i + 1] << std::endl;
                    ip.h_mesh = atof(argv[++i]);
                    continue;
                }
                if (strcmp(argv[i], "-R0") == 0 || strcmp(argv[i], "-R0=") == 0){
                    std::cout << "initial radius = " << argv[i + 1] << std::endl;
                    ip.R0 = atof(argv[++i]);
                    continue;
                }
                if (strcmp(argv[i], "-m") == 0 || strcmp(argv[i], "--model") == 0){
                    std::cout << "model number = " << argv[i + 1] << std::endl;
                    ip.ModelNum = atoi(argv[++i]);
                    continue;
                }
            }
            return ip;
        }
    };
    InputParams params = InputParams::processMainArgs(argc, argv);
    std::vector<int> vrtmaterial;
    net_t mesh = generate_eigth_sphere_part(params.R0, params.h_mesh, vrtmaterial);
    double thickness = 0.5; //mm
    set_thickness(mesh, thickness);

    net_t_set_elems_neighbours(mesh);		//order is unimportant
    nets_t _nets = nets_t_get_net(1); _nets.nets[0] = mesh;
    set_relax_consts(_nets, 0.02, 0.05);		//0.02, 0.05
    nets_t_set_relax_state(_nets, point_t_get_point(1, 0, 0));
    double P = 80;  // Hg mm
    double delta = 1e-7;
    double err = 1e-12;
    double maxits = 100000;
    to_stl(_nets, "../bin/preresult");
    //exit(0);
    set_default_elastic_model(EMOD_NEOGOOK);
    double R_exact =predicted_radius(P, params.R0, thickness);
    test_compute_nets(_nets, P, delta, maxits, err, 100, vrtmaterial);
    to_stl(_nets, "../bin/result");
    std::cout << "Exact radius = " << R_exact << std::endl;
    double R_model = test_spheriable(mesh);
    std::cout << "|R_exact - R_eval| = " << fabs(R_exact - R_model) << std::endl;

    return 0;
}

double predicted_radius(double P /*Hg mm*/, double R0 /*mm*/, double T /*Pa*/){
    double mu = 1.0e6 / 3;
    P *= 133.3;
    double k = P * R0 / (2 * mu * T);
    const int n = 8;
    double a[n] = { -k, 1, 0, 0, 0, 0, 0, -1 };
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
        if ((0 < z[2*i] && z[2*i] < 1)) {
            real = i;
            rreal = z[2*i + 1];
            if (fabs(rreal) < eps)
                solved = 1;
            break;
        }

    for (; i < n - 1; ++i)
        if ((0 < z[2*i] && z[2*i] < 1) && (fabs(z[2*i+1]) < fabs(rreal)) && (fabs(z[2*i+1]) < eps)){
            real = i;
            rreal = z[2*i + 1];
            solved = 2;
        }
    if (!solved || fabs(z[real * 2 + 1]) > eps)
        throw std::runtime_error("Failed to solve the equation");

    return R0 / z[2*real];
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
    std::cout << "Evaluate radius = " << r << std::endl;
    std::cout << "Maximal deviation from radius is " << max_div << std::endl;
    std::cout << "Middle deviation from radius is " << mid_div << std::endl;

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

void test_compute_nets(nets_t nets, double P, double delta, unsigned int n, double eps, int freq, std::vector<int>& vrt_mask){
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
                            printf("current possition: "); point_t_dump(cam->loc); printf("\n");
                            printf("camPitch = %lg, camYaw = %lg\n", cam->camPitch, cam->camYaw);
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
                for (int p = 0; p < 3; ++p) {
                    if (!mask) break;
                    if (mask % 2) nets.nets[0].vrtx.nodes[j]->next.coord[p] = 0;
                    mask /= 2;
                }
            }

            if (!flag) {
                double diviation = get_diviation(nets);
                printf("it = %d, diviation = %lg, relation = %lg\n", i, diviation / delta, diviation / init_div);
            }

            update_Mid_Shift(nets);

            if (i == 1) {
                double max_Mid_Shift = 0;
                unsigned int l;
                for (l = 0; l < get_Cnt_Nodes(); l++) {
                    double len = SQR_LEN(get_Mid_Shift(l));
                    if (max_Mid_Shift < len) max_Mid_Shift = len;
                }
                printf("Max_Shift = %e\n", sqrt(max_Mid_Shift));
                init_div = get_diviation(nets);
            }
            if (!(i % freq)) {
                res = get_Mid_Shift_len(freq);
                double rms = res / get_Cnt_Nodes(); // delta ;
                printf("RMS shift per iter of node = %e mm, relation = %lg, cr = %d\n", rms, res / init_div, crush);
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
    printf("Made %u iterations, relation = %lg\n", --i, res / init_div);

}

void set_thickness(net_t& mesh, double thickness){
    for (int i = 0; i < mesh.vrtx.count; ++i)
        mesh.vrtx.nodes[i]->h = thickness;
}


//mesh generator functions//
static double global_radius = 1.0;
static double global_size = 0.1;

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
        vertexmaterial[edge[2*i+0]-1] = edgematerial[i];
        vertexmaterial[edge[2*i+1]-1] = edgematerial[i];
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

net_t generate_eigth_sphere_part(double R0, double size, std::vector<int>& vm){
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

    printf("size = %lf\n", size);

    printf("\n * Generating surface mesh\n"); /* just an example */
    create_sphere_quarter(&nV, vertex, &nF, face, facematerial, &nE, edge, edgematerial,  nnV, nnF, nnE,  R0, size);
    mark_vertices(nV, vertexmaterial, nE, edge, edgematerial, 0);
    // Write surface mesh
    write_front_gmv   ("../bin/surface.gmv", nV, vertex, nF, face, facematerial);
    write_front       ("../bin/surface.smv", nV, vertex, nF, face, facematerial); // for smv
    write_custom_format("../bin/surface.txt", nV, vertex, vertexmaterial, nF, face);

    printf("INFO: nV = %d, nF = %d, nE = %d\n", nV, nF, nE);

    /* Add layers */
    printf("\n * Adding layers\n");

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