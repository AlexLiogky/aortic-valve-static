#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <math.h>
#include <string.h>
#include <sys/time.h>

#include <SDL/SDL.h>
#include <GL/gl.h>
#include <GL/glu.h>

#include "nets.h"
#include "stress.h"
#include "precomputation.h"
#include "computation.h"
#include "postcomputation.h"
#include "intersection.h"
#include "save-data.h"
#include "format_in.h"
#include "separate.h"
#include "sew_leaf.h"
#include "bound-box.h"
#include "camera.h"
//#include "Solver.h"
#include "World.h"

#define RES_STL         "results/sew0"
#define RES_LEAF_STL    "results/sew0_leaf"
#define RES_A_L_STL     "results/sew0_aorta_leaf"
#define INPUT           "data/mesh-templates/templ-17-800.txt"//"templ-25-100.txt" "leaflets/leaf-1200.txt" "separate2" "leaflets.txt"
#define OUTPUT          "results/full-res4"
#define S_STL           "results/result11"
#define LEAF_STL        "results/res0"
#define AORTA_IN        "data/aorta/aorta.nts"
#define BND_IN          "data/aorta/aorta_rec_bnd.bnd"

#include <string>
#include <sstream>
#include <fstream>
#include <iostream>
using namespace std;


std::ofstream gLog;

//######################################################################
long get_msec_time(struct timeval start, struct timeval end){
	long seconds  = end.tv_sec  - start.tv_sec;
    long useconds = end.tv_usec - start.tv_usec;
	long mtime = ((seconds) * 1000 + useconds / 1000.0) + 0.5;
	return mtime;
}

long double get_msec_time_(struct timeval start, struct timeval end){
	long double seconds  = end.tv_sec  - start.tv_sec;
    long double useconds = end.tv_usec - start.tv_usec;
	long double mtime = (seconds) * 1000 + useconds / 1000.0;
	return mtime;
}

void print_max_min(nets_t nets){
	int nets_c = nets.count, i;
	double max = 0, min = nets.nets[0].springs.springs[0][0].l_0;
	for(i = 0; i < nets_c; i++){
		net_t net = nets.nets[i];
		int s_c = net.springs.count, j;
		for (j = 0; j < s_c; j++){
			double len = net.springs.springs[j][0].l_0;
			if (max < len) max = len;
			if (min > len) min = len;
		}
	}
	printf("max-spr = %lg, min-spr = %lg\n", max, min);
}

void save_data(nets_t nets, char* file_name){
	char name[1024] = {0};
	strcpy(name, file_name);
	int len = strlen(file_name);
	nets_t net = nets_t_get_net(1);
	for (unsigned char i = 0; i < nets.count; i++){
		name[len] = '0' + i;
		net.nets[0] = nets.nets[(int)i];
		printf("name = <%s>\n", name);
		to_stl(net, name);
	}
}

nets_t download_aorta(char* file_name){
	nets_t aorta = download_nets_from_file(file_name);
	precomputation(aorta);
	return aorta;
}

nets_t get_system(nets_t aorta, char* bnd, point_t blood_direction, char* leaf1, char* leaf2, char* leaf3){
	point_t shift = blood_direction;

	point_t init_points[2], final_points[2];
	point_t init_shift = point_t_get_point(0, -9, 0);
	point_t fiber_direction = point_t_get_point(0, 1, 0);

	point_t commissur[3];
	bnds_t bnds = read_bnds(bnd);
	for (int i = 0; i < bnds.cnt; i++)
		commissur[i] = bnds.bnds[i].line[0];
	/*point_t bottom[3] = {{{-11.941, -36.7436, 2.81578}}, {{-3.33692, -17.0099, 14.5751}}, {{-14.1942, -21.4188, -4.19314}}};
	sew_line_t* sews = init_sew_line_on_aorta(aorta, commissur, bottom, 3);
	//sew_line_t_plate_dump(&sews[0]);
	point_t p = {{-7.4382, -28.391, 20.944}};
	sew_line_t_add_point(p, &sews[1]);
	sew_line_t_add_point(VEC(-4.994, -26.491, 18.874), &sews[1]);
	sew_line_t_add_point(VEC(-3.515, -23.861, 18.608), &sews[1]);
	sew_line_t_add_point(VEC(-16.265, -13.279, 10.909), &sews[1]);
	bnds_t newbnds = {3, calloc(3, sizeof(line_t))};
	for (int i = 0; i < 3; i++)
		newbnds.bnds[i] = sew_line_t_get_line_on_aorta(&sews[i]);
	printf("count points in leaf[0] = %d\n", newbnds.bnds[1].pnt_cnt);
	//newbnds.bnds[0] = sew_line_t_get_line_on_aorta(&sews[0]);
	line_t_dump(newbnds.bnds[1]);
	bnds = newbnds;*/


	nets_t leaflet1 = formated_in(leaf2);
	nets_t_set_relax_state(leaflet1, fiber_direction);
	to_stl(leaflet1, (char*)LEAF_STL);
	init_points[0] = leaflet1.nets[0].vrtx.nodes[1][0].coord;
	init_points[1] = leaflet1.nets[0].vrtx.nodes[4][0].coord;
	final_points[0] = commissur[1];
	final_points[1] = commissur[0];
	sew_leaflet_to_aorta(leaflet1.nets[0], init_points, init_shift,\
						aorta.nets[0], final_points, shift, 4, bnds.bnds[0]);

	nets_t leaflet2 = formated_in(leaf1);
	nets_t_set_relax_state(leaflet2, fiber_direction);
	init_points[0] = leaflet2.nets[0].vrtx.nodes[1][0].coord;
	init_points[1] = leaflet2.nets[0].vrtx.nodes[4][0].coord;
	final_points[0] = commissur[2];
	final_points[1] = commissur[1];
	sew_leaflet_to_aorta(leaflet2.nets[0], init_points, init_shift,\
						aorta.nets[0], final_points, shift, 4, bnds.bnds[1]);

	nets_t leaflet3 = formated_in(leaf3);
	nets_t_set_relax_state(leaflet3, fiber_direction);
	init_points[0] = leaflet3.nets[0].vrtx.nodes[1][0].coord;
	init_points[1] = leaflet3.nets[0].vrtx.nodes[4][0].coord;
	final_points[0] = commissur[0];
	final_points[1] = commissur[2];
	sew_leaflet_to_aorta(leaflet3.nets[0], init_points, init_shift,\
						aorta.nets[0], final_points, shift, 4, bnds.bnds[2]);

	int n_objects = 4;
	nets_t leaflet = nets_t_get_net(n_objects);
	leaflet.nets[0] = leaflet1.nets[0];
	leaflet.nets[1] = leaflet2.nets[0];
	leaflet.nets[2] = leaflet3.nets[0];
	leaflet.nets[3] = aorta.nets[0];
	net_t_set_state(&leaflet.nets[3], 1);
	precomputation(leaflet);

	return leaflet;
	//return aorta;
}

nets_t get_valve_from_system(nets_t system){
	nets_t valve = nets_t_get_net(system.count - 1);
	for (unsigned int i = 0; i < system.count - 1; i++)
		valve.nets[i] = system.nets[i];

	return valve;
}

void print_nets_statistic(nets_t leaflet, point_t shift){
                                        set_contact_recognition_resolution(1.1 * 0.1);
	printf("Contact_Resolution = %lg\n", get_Contact_Resolution());
        gLog << "Contact Resolution = " << get_Contact_Resolution()<< endl;
	nets_t curleaflet = leaflet;//get_valve_from_system(leaflet);
	print_statistic(curleaflet);
	double l_free[10];
	for (unsigned int i = 0; i < leaflet.count; ++i)
        l_free[i] = net_t_get_len_free_edge(leaflet.nets[i]);

	double h = nets_t_get_coapt_depth(curleaflet, &shift);
	int _pow = -1;
	                                    set_contact_recognition_resolution(1.1 * 0.35);
	double h_c = nets_t_get_coapt_intersect_depth(curleaflet, &_pow);
                                        set_contact_recognition_resolution(1.1 * 0.1);
    printf("l_free[leaf] = {%lg, %lg, %lg}\n", l_free[0], l_free[1], l_free[2]);
	printf("%s: h = %lg, h_c = %lg\n", OUTPUT, h, h_c);
        gLog << "l_free[leaf] = { " <<  l_free[0] << ", " << l_free[1] << ", " << l_free[2] << " }\n";
        gLog << "0-th level:\n" << " h = " << h << "\n";
        gLog << " h_c = " << h_c << ", count of central points = " << _pow << endl;

	for (int i = 0; i < 1; i++){
		curleaflet = create_next_hierarchical_nets(curleaflet);
            gLog << i + 1 << "-th level:\n";
		//auto_set_contact_recognition_consts(curleaflet, NULL);

		                                    set_contact_recognition_resolution(1.1 * 0.35);
		double h_c = nets_t_get_coapt_intersect_depth(curleaflet, &_pow);
		                                    set_contact_recognition_resolution(1.1 * 0.1);
		double h = nets_t_get_coapt_depth(curleaflet, &shift);
            gLog << " h = " << h << "\n";
            gLog << " h_c = " << h_c << ", count of central points = " << _pow << endl;
		double f_area[3] = {};
		double S_max = 0, S_min = 0;
		printf("Compute full leaflet coaptation area:\n");
            gLog << " Full leaflet coaptation area:\n";
		for (int i = 0; i < 3; i++){
			 f_area[i] = net_t_get_coapt_area(curleaflet.nets[i]);
			 if (i == 0) S_max = f_area[i], S_min = f_area[i];
			 S_max = (S_max < f_area[i]) ? f_area[i] : S_max;
			 S_min = (S_min > f_area[i]) ? f_area[i] : S_min;
			 printf("  full area[%d] = %lg\n", i, f_area[i]);
                gLog << "  full area[" << i << "] = " << f_area[i] << endl;
		 }
		double p_area[6] = {};
		printf("Compute partition leaflet coaptation area:\n");
            gLog << " Partition leaflet coaptation area:\n";
		for (int i = 0; i < 3; i++){
			p_area[2 * i] = nets_t_get_coapt_area_from_to(curleaflet, i%3, (i+1)%3);
			p_area[2 * i + 1] = nets_t_get_coapt_area_from_to(curleaflet, (i+1)%3, i%3);
			printf("partion area[%d -> %d] = %lg\n", i%3, (i+1)%3, p_area[2 * i]);
			printf("partion area[%d -> %d] = %lg\n", (i + 1)%3, (i)%3, p_area[2 * i + 1]);
                gLog << "  partion area[" << i%3 << " -> " << (i+1)%3 << "] = " << p_area[2 * i] << "\n";
                gLog << "  partion area[" << (i + 1)%3 << " -> " << (i)%3 << "] = " << p_area[2 * i + 1] << endl;
		}
        double h_mid[10];
        for (unsigned int i = 0; i < leaflet.count; ++i)
            h_mid[i] = (p_area[2 * ((i + 2) % 3) + 0] + p_area[2 * i + 1]) / l_free[i];
                gLog << " h_mid[leaf] = { " <<  h_mid[0] << ", " << h_mid[1] << ", " << h_mid[2] << " }\n";
		printf("%s: h = %lg, h_c = %lg, S_max = %lg, S_min = %lg\n", OUTPUT, h, h_c, S_max, S_min);
		printf("%s: h_mid[leaf] = { %lg, %lg, %lg }\n", OUTPUT, h_mid[0], h_mid[1], h_mid[2]);
	}
}

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
    const double max_deform = 0.02;
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


void display(world_t* world, camera_t* cam)
{
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
	glLoadIdentity();
	camera_t_Control(cam);
	//drawSkybox(50);
	camera_t_UpdateCamera(cam);
	for (unsigned i = 0; i < world->dynamic_nets.count; ++i)
        renderSoftBody(world->dynamic_nets.nets[i], point_t_get_point(0.296, 0.221, 0.231));

    for (unsigned i = 0; i < world->static_nets.count; ++i)
        renderSoftBody(world->static_nets.nets[i], point_t_get_point(205.0/256, 200.0/256, 239.0/256));

}

template <class TT>
void displayT(TT* s, camera_t* cam)
{
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
	glLoadIdentity();
	camera_t_Control(cam);
	//drawSkybox(50);
	camera_t_UpdateCamera(cam);
	nets_t& nets = s->getDynamicNets();
	for (unsigned i = 0; i < nets.count; ++i)
        renderSoftBody(nets.nets[i], point_t_get_point(0.296, 0.221, 0.231));

    nets_t& nets1 = s->getStaticNets();
    for (unsigned i = 0; i < nets1.count; ++i)
        renderSoftBody(nets1.nets[i], point_t_get_point(205.0/256, 200.0/256, 239.0/256));

}

void run_comupation(world_t* world){
    SDL_Init(SDL_INIT_EVERYTHING);
	SDL_SetVideoMode(1280,800,32,SDL_OPENGL);
	Uint32 _start;
	SDL_Event event;
	int running=1;
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

	while(running)
	{
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
		compute_nets_time(1000.0/60, world, 1000);
        display(world, cam);
		SDL_GL_SwapBuffers();
		if(1000.0/60>SDL_GetTicks()-_start)
			SDL_Delay(1000.0/60-(SDL_GetTicks()-_start));
	}

	SDL_Quit();
	camera_t_destruct(cam);
}

template<class TT>
void run_comupationT(TT* s, unsigned int maxtime_sec = 3600 * 24){
    SDL_Init(SDL_INIT_EVERYTHING);
	SDL_SetVideoMode(1280,800,32,SDL_OPENGL);
	Uint32 _start, start;
	SDL_Event event;
	int running=1;
	int compute=1;
	int render_mode=1;
	float angle=45;
	glClearColor(0,0,0,1);
	glMatrixMode(GL_PROJECTION);
		glLoadIdentity();
		gluPerspective(angle,1280.0/800.0,1,1000);
	glMatrixMode(GL_MODELVIEW);
	glEnable(GL_DEPTH_TEST);
	camera_t* cam = camera_t_construct(point_t_get_point(16.0528, 27.7954, -83.2535));//-45.8, -22, 34.6));
	cam->camPitch = -13;//-8;
	cam->camYaw = 174;//308;
	//point_t_dump(camera_t_getVector(cam));

	start = SDL_GetTicks();
	unsigned int maxtime_msec = maxtime_sec * 1000;
	while(running && SDL_GetTicks() - start < maxtime_msec)
	{
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
                            compute ^= 1;
							break;
                        }
                        case SDLK_EQUALS:
                        {
                            solver_t data = s->getSolverData();
                            data.delta /= 0.9;
                            s->setSolverData(data);
                            printf("delta = %lg\n", data.delta);
                                gLog << "new delta = " << data.delta << "\n";
                            break;
                        }
                        case SDLK_MINUS:
                        {
                            solver_t data = s->getSolverData();
                            data.delta *= 0.9;
                            s->setSolverData(data);
                            printf("delta = %lg\n", data.delta);
                                gLog << "new delta = " << data.delta << "\n";
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
		double dt = 1000.0/60;
		if (compute)
            s->compute_nets_time(dt, 1000);
        displayT<TT>(s, cam);
		SDL_GL_SwapBuffers();
		if(dt>SDL_GetTicks()-_start)
			SDL_Delay(dt-(SDL_GetTicks()-_start));
	}

	SDL_Quit();
	camera_t_destruct(cam);
}


//######################################################################
//######################################################################
int main(int argc, char* argv[]){
    std::string directory = "/home/alex/Desktop/MV/valve model-static/Aortic-Valve/src/additional/Haj-mesh/template-20/";
    std::string filename = "haj-valve";
    std::string res_fname = filename + "-res";
    std::string res_dir = "haj-results/";
    int nres = 1;
    long int tm = time(NULL);
    res_fname += "-" + to_string(tm);
    std::string logfile = filename + "-" + to_string(tm) + "_log.txt";

    switch (argc)
    {
        case 1: break;
        case 6: logfile = string(argv[5]);
        case 5: res_dir = string(argv[4]);
                if (res_dir[res_dir.length() - 1] != '/')
                    res_dir += "/";
        case 4: res_fname = string(argv[3]);
        case 3: filename = string(argv[2]);
        case 2: directory = string(argv[1]);
                if (directory[directory.length() - 1] != '/')
                    directory += "/";
                break;
        default:
            cout << "Wrong input\n";
            cout << "Usage: " << argv[0] << " <directory_name> <fname> <result file name> <result directory> <logfile name>\n\n";
            return 0;
    }

    switch (argc)
    {
        case 1: cout << "Input directory: " << directory << "\n";
        case 2: cout << "Input file name: " << filename << "\n";
        case 3: cout << "Result file name: " << res_fname + "-<#res>\n";
        case 4: cout << "Result directory: " << res_dir << "\n";
        case 5: cout << "Logfile: " << res_dir + logfile << "\n";
        default: break;
    }

    gLog.open(res_dir + logfile, ios::trunc);

        gLog << "Open " << directory << filename << endl;
    nets_t test = formated_in((directory + filename).c_str());//"/home/alex/Desktop/MV/valve model-static/Aortic-Valve/src/additional/Labrosse-mesh-generator/sample");
    for (unsigned int i = 0; i < test.count; ++i)
        for (unsigned int j = 0; j < test.nets[i].vrtx.count; ++j)
            test.nets[i].vrtx.nodes[j]->h = 0.3;

        gLog << "Elems per net = " << test.nets[0].elems.count << "\n";
        gLog << "Nodes per net = " << test.nets[0].vrtx.count << "\n";
        gLog << "Springs per net = " << test.nets[0].springs.count << "\n";
        gLog << "Initial full area per net = " << net_t_get_full_area(test.nets[0]) << "\n";
    nets_t static_nets = nets_t_get_net(0);
    nets_t dynamic_nets = test;
        //to_stl(dynamic_nets, result.c_str()); return 0;

    nets_t_set_relax_state(test, point_t_get_point(0, 1, 0));
    precomputation(test);
    printf("len = %lg\n", net_t_get_len_free_edge(test.nets[0]));
        gLog << "Initial free edge len per net = " << net_t_get_len_free_edge(test.nets[0]) << "\n";
        gLog << "Fixed len net = " << net_t_get_len_fix_edge(test.nets[0]) << endl;
    solver_t solver_data = {25e-7, 0.001};
        gLog << "delta = " << solver_data.delta << "\n";
    wrld_cnd_t conditions = {80.0};
        gLog << "P = " << conditions.P << endl;
    World s(dynamic_nets, static_nets, conditions, solver_data, 200*Allow_shift, 200*Max_shift, 0.05);
    update_nets(test);
    struct timeval start1, end1;
    gettimeofday(&start1, NULL);
    run_comupationT<World>(&s, 3600 * 2);
    gettimeofday(&end1, NULL);
    printf("Time of computation = %ld ms\n", get_msec_time(start1, end1));
        gLog << "Final full area per net = " << net_t_get_full_area(test.nets[0]) << "\n";
        gLog << "Final free edge len per net = " << net_t_get_len_free_edge(test.nets[0]) << endl;
        gLog << "Time of computation = " << get_msec_time(start1, end1) << endl;

    to_stl(dynamic_nets, (res_dir + res_fname + "-" + to_string(nres++)).c_str());

    nets_t dynamic_nets1 = create_next_hierarchical_nets(dynamic_nets);
    precomputation(dynamic_nets1);
    World s1(dynamic_nets1, static_nets, conditions, solver_data, 200*Allow_shift, 200*Max_shift, 0.03);
    update_nets(test);
    gettimeofday(&start1, NULL);
    run_comupationT<World>(&s1, 3600);
    gettimeofday(&end1, NULL);
    printf("Time of computation = %ld ms\n", get_msec_time(start1, end1));
        gLog << "Final full area per net = " << net_t_get_full_area(dynamic_nets1.nets[0]) << "\n";
        gLog << "Final free edge len per net = " << net_t_get_len_free_edge(dynamic_nets1.nets[0]) << endl;
        gLog << "Time of computation = " << get_msec_time(start1, end1) << endl;

    to_stl(dynamic_nets1, (res_dir + res_fname + "-" + to_string(nres++)).c_str());
    print_nets_statistic(dynamic_nets1, point_t_get_point(0, 0, 1));

//    collision_t collision_data = collision_data_t_construct(dynamic_nets, static_nets, 10);
//    wrld_cnd_t conditions = {80.0};
//    world_t* world = world_t_construct(dynamic_nets, static_nets, solver_data , collision_data, conditions);
//    set_initial_solving_params(world);
//    struct timeval start1, end1;
//	gettimeofday(&start1, NULL);
//    run_comupation(world);
//    gettimeofday(&end1, NULL);
//	printf("Time of computation = %ld ms\n", get_msec_time(start1, end1));

    return 0;
//	nets_t aorta = download_aorta((char*)AORTA_IN);
//	printf("elems.count = %u\n", aorta.nets[0].elems.count);
//	printf("Contact_Resolution = %lg\n", get_Contact_Resolution());
//	point_t shift = point_t_get_point(5, 1, -3.3);
//	nets_t leaflet = get_system(aorta, (char*)BND_IN, shift, (char*)INPUT, (char*)INPUT, (char*)INPUT);
//
//    /*nets_t dynamic_nets = nets_t_get_net(3);
//    for (int i = 0; i < 3; ++i) dynamic_nets.nets[i] = leaflet.nets[i];
//    nets_t static_nets = nets_t_get_net(1);
//    static_nets.nets[0] = leaflet.nets[3];
//    solver_t solver_data = {2e-7, 0.001};
//    collision_t collision_data = collision_data_t_construct(dynamic_nets, static_nets, 10);
//    wrld_cnd_t conditions = {80.0};
//    world_t* world = world_t_construct(dynamic_nets, static_nets, solver_data , collision_data, conditions);
//    set_initial_solving_params(world);
//    struct timeval start1, end1;
//	gettimeofday(&start1, NULL);
//    run_comupation(world);
//    gettimeofday(&end1, NULL);
//	printf("Time of computation = %ld ms\n", get_msec_time(start1, end1));*/
//
//    nets_t dynamic_nets = nets_t_get_net(3);
//    for (unsigned int i = 0; i < dynamic_nets.count; ++i) dynamic_nets.nets[i] = leaflet.nets[i];
//    //dynamic_nets.nets[0] = leaflet.nets[1];
//    nets_t static_nets = nets_t_get_net(1);
//    static_nets.nets[0] = leaflet.nets[3];
//    solver_t solver_data = {15e-7, 0.001};
//    wrld_cnd_t conditions = {80.0};
//    World s(dynamic_nets, static_nets, conditions, solver_data, 200*Allow_shift, 200*Max_shift);
//    struct timeval start1, end1;
//	gettimeofday(&start1, NULL);
//    run_comupationT<World>(&s);
//    gettimeofday(&end1, NULL);
//	printf("Time of computation = %ld ms\n", get_msec_time(start1, end1));
//	to_stl(leaflet, (char*)OUTPUT);
//    return 0;
//
//	struct timeval start, end;
//	gettimeofday(&start, NULL);
//	double P = 80; //mm Hg
//	double delta = 2e-7;//3e-7;
//	printf("delta = %lg\n", delta);
//	printf("Contact_Resolution = %lg\n", get_Contact_Resolution());
//	int max_nsteps = 28000, freq = 150;
//	double eps = 0;
//	//compute_nets(leaflet, P, delta, max_nsteps, eps, freq);
//	nets_t curleaflet = get_valve_from_system(leaflet);
//	curleaflet = create_next_hierarchical_nets(curleaflet);
//	int jj = 3;
//	for (unsigned int ii = 0; ii < leaflet.nets[jj].elems.count; ++ii){
//		leaflet.nets[jj].elems.elems[ii]->coef *= -1;
//	}
//	to_stl(leaflet, (char*)OUTPUT);
//	gettimeofday(&end, NULL);
//	//printf("Time of precomputation = %ld ms\n", get_msec_time(start, end));
//
//	print_nets_statistic(leaflet, shift);
//
//	printf("P = %lg, %s\n", P, OUTPUT);
//
//	//world_t_destruct(world);
//	//collision_data_t_destruct(collision_data);
//	nets_t_destruct(curleaflet);
//	nets_t_surfacial_free(dynamic_nets);
//	nets_t_surfacial_free(static_nets);
//	nets_t_surfacial_free(aorta);
//	nets_t_destruct(leaflet);
//
//
//	return 0;
}
//######################################################################
//######################################################################
/*TODO:
 * т.к. теперь есть exact_box, то можно создать замену для функции box_t_local_node_to_net_projection(...)
 * сделать проверку правильности последовательности (по state) в compute_nets(...)
 */

 /* можно отключить update_node у неподвижных вершин */
 /*центральная кооптация считается неверно, нужно понять, почему*/

 /*изменена контактная сила*/

/*
 * создать считыватель всех параметров из config.conf
 * сделать полный рефакторинг кода
 * */

/*
 * имеет смысл делать константы сил контакта
 * пропорциональными sqrt(S) ~ contact_resolution
 * k = k * get_contact_resolution()
 */

/*nets_t leaflet = formated_in(INPUT);
	nets_t_set_relax_state(leaflet, point_t_get_point(0, 0, 1));
	point_t shift = point_t_get_point(0, 0, 1);

	double P = 1; //mm Hg
	double delta = 3e-7;//3e-7;
	printf("delta = %lg\n", delta);
	int max_nsteps = 3000, freq = 150;
	double eps = 0;
	//compute_nets(leaflet, P, delta, max_nsteps, eps, freq);
	precomputation(leaflet);
	P = 80, max_nsteps = 24000;
	compute_nets(leaflet, P, delta, max_nsteps, eps, freq);
	leaflet = create_next_hierarchical_nets(leaflet);
	leaflet = create_next_hierarchical_nets(leaflet);
	leaflet = create_next_hierarchical_nets(leaflet);

	printf("%s: coapt_depth = %lg\n", RES_STL, nets_t_get_coapt_depth(leaflet, &shift));

	to_stl(leaflet, "result");*/

