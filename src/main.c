#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <math.h>
#include <string.h>
#include <sys/time.h>

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

#define RES_STL         "results/sew0"
#define RES_LEAF_STL    "results/sew0_leaf"
#define RES_A_L_STL     "results/sew0_aorta_leaf"
#define INPUT           "data/mesh-templates/templ-17-800.txt"//"templ-25-100.txt" "leaflets/leaf-1200.txt" "separate2" "leaflets.txt"
#define OUTPUT          "results/full-res4"
#define S_STL           "results/result11"
#define LEAF_STL        "results/res0"
#define AORTA_IN        "data/aorta/aorta.nts"
#define BND_IN          "data/aorta/aorta_rec_bnd.bnd"


//######################################################################
long get_msec_time(struct timeval start, struct timeval end){
	long seconds  = end.tv_sec  - start.tv_sec;
    long useconds = end.tv_usec - start.tv_usec;
	long mtime = ((seconds) * 1000 + useconds / 1000.0) + 0.5;
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
	printf("Contact_Resolution = %lg\n", get_Contact_Resolution());
	nets_t curleaflet = get_valve_from_system(leaflet);
	print_statistic(curleaflet);
	double h = nets_t_get_coapt_depth(curleaflet, &shift);
	double h_c = nets_t_get_coapt_intersect_depth(curleaflet);
	printf("%s: h = %lg, h_c = %lg\n", OUTPUT, h, h_c);

	for (int i = 0; i < 2; i++){
		curleaflet = create_next_hierarchical_nets(curleaflet);
		//auto_set_contact_recognition_consts(curleaflet, NULL);
		double h_c = nets_t_get_coapt_intersect_depth(curleaflet);
		double h = nets_t_get_coapt_depth(curleaflet, &shift);
		double f_area[3] = {};
		double S_max = 0, S_min = 0;
		printf("Compute full leaflet coaptation area:\n");
		for (int i = 0; i < 3; i++){
			 f_area[i] = net_t_get_coapt_area(curleaflet.nets[i]);
			 if (i == 0) S_max = f_area[i], S_min = f_area[i];
			 S_max = (S_max < f_area[i]) ? f_area[i] : S_max;
			 S_min = (S_min > f_area[i]) ? f_area[i] : S_min;
			 printf("full area[%d] = %lg\n", i, f_area[i]);
		 }
		double p_area[6] = {};
		printf("Compute partition leaflet coaptation area:\n");
		for (int i = 0; i < 3; i++){
			p_area[2 * i] = nets_t_get_coapt_area_from_to(curleaflet, i%3, (i+1)%3);
			p_area[2 * i + 1] = nets_t_get_coapt_area_from_to(curleaflet, (i+1)%3, i%3);
			printf("partion area[%d -> %d] = %lg\n", i%3, (i+1)%3, p_area[2 * i]);
			printf("partion area[%d -> %d] = %lg\n", (i + 1)%3, (i)%3, p_area[2 * i + 1]);
		}
		printf("%s: h = %lg, h_c = %lg, S_max = %lg, S_min = %lg\n", OUTPUT, h, h_c, S_max, S_min);
	}
}
//######################################################################
//######################################################################
int main(){
	nets_t aorta = download_aorta((char*)AORTA_IN);
	//to_stl(aorta, "res-aorta");
	printf("elems.count = %u\n", aorta.nets[0].elems.count);
	//net_t_set_state(&aorta.nets[0], 1);
	printf("Contact_Resolution = %lg\n", get_Contact_Resolution());
	//set_contact_recognition_resolution(0.51098);
	/*box_t box = box_t_construct(aorta, get_Contact_Resolution());
	nets_t_update_box(aorta, &box);
	point_t intersect = {};
	point_t line[2];
	line[0] = VEC(-32.2005, -36.3946, 8.15829);
	//printf("find "); point_t_dump(line[0]);
	line[1] = SUM(line[0], VEC(1, 1, 1));
	if (line_to_boxed_nets_intersection(line, aorta, 0, box, &intersect)){
		printf("Intersection at ");
		point_t_dump(intersect);
	}*/
	point_t shift = point_t_get_point(5, 1, -3.3);
	nets_t leaflet = get_system(aorta, (char*)BND_IN, shift, (char*)INPUT, (char*)INPUT, (char*)INPUT);
	/*nets_t researched = nets_t_get_net(2);
	nets_t leaf = nets_t_get_net(1);
	researched.nets[0] = leaflet.nets[1];
	leaf.nets[0] = leaflet.nets[1];
	to_stl(leaf, RES_LEAF_STL);
	researched.nets[1] = leaflet.nets[3];
	to_stl(researched, RES_A_L_STL);*/
	//to_stl(leaflet, RES_STL);

	struct timeval start, end;
	gettimeofday(&start, NULL);
	double P = 80; //mm Hg
	double delta = 2e-7;//3e-7;
	printf("delta = %lg\n", delta);
	printf("Contact_Resolution = %lg\n", get_Contact_Resolution());
	int max_nsteps = 28000, freq = 150;
	double eps = 0;
	compute_nets(leaflet, P, delta, max_nsteps, eps, freq);
	nets_t curleaflet = get_valve_from_system(leaflet);
	curleaflet = create_next_hierarchical_nets(curleaflet);
	/*for (int i = 0; i < 3; i++) leaflet.nets[i] = curleaflet.nets[i];
	auto_set_contact_recognition_consts(curleaflet, NULL);
	max_nsteps = 4000;
	for (int i = 1; i < 2; i++){
		compute_nets(leaflet, P, delta, max_nsteps, eps, freq);
		printf ("MADE %d iterations\n", (i + 1) * max_nsteps);
		char res[10] = {};
		sprintf(res, "res%d", i+1);
		//to_stl(curleaflet, res);
	}*/
	//for (int jj = 0; jj < 4; ++jj)
	int jj = 3;
	for (unsigned int ii = 0; ii < leaflet.nets[jj].elems.count; ++ii){
		leaflet.nets[jj].elems.elems[ii]->coef *= -1;
	}
	to_stl(leaflet, (char*)OUTPUT);
	gettimeofday(&end, NULL);
	printf("Time of precomputation = %ld ms\n", get_msec_time(start, end));

	print_nets_statistic(leaflet, shift);

	printf("P = %lg, %s\n", P, OUTPUT);

	return 0;
}
//######################################################################
//######################################################################
/*TODO:
 * т.к. теперь есть exact_box, то можно создать замену для функции box_t_local_node_to_net_projection(...)
 * сделать проверку правильности последовательности (по state) в compute_nets(...)
 */

 /* можно отключить update_node у неподвижных вершин */

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

