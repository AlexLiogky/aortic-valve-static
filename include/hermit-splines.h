#ifndef _HSPLINES_H
#define _HSPLINES_H

#include "point.h"

#ifdef __cplusplus
extern "C" {
#endif

typedef struct p2d{
	double pnt[2];
} p2d;

typedef struct spdot_t{
	p2d pnt;	//point
	p2d der;	//derivation
	struct spdot_t* next;
	struct spdot_t* prev;
}spdot_t;

typedef struct spline_t{
	int len;
	spdot_t* head;
}spline_t;


void p2d_dump(p2d obj);				//debug print of p2d
void spdot_t_dump(spdot_t* obj);		//debug print of spdot_t
void spline_t_full_dump(spline_t sp); 	//debug print of spline_t
spline_t spline_construct(spdot_t* head); 	//constructor
void spdot_destruct(spdot_t* dot); 			//destructor
p2d get_p2d(double pnt[2]); 				//constructor
spdot_t* spdot_t_construct(double pnt[2], double der[2]); //constructor
void spdot_t_set_der(spdot_t* dot, double der[2]); //set derivation "der" into "dot"
void spline_destruct(spline_t spline); 		//destructor
void spline_t_add_dot(spline_t* spline, spdot_t* insert, spdot_t* prev); //insert "insert" after "prev" into "spline"
point_t extend_2d_to_3d_point(p2d p);
p2d cut_3d_to_2d_point(point_t p);
double sqr_dist_to_2D_line_fragment(p2d frag[2], p2d p); //compute distance from "p" to line fragment with ends -"frag"

//return pointer on the nearest dot in spline to "pnt", if "id" set id = returned spdot_id
spdot_t* get_nearest_prev_dot(spline_t spline, double pnt[2], int* id);

//automatically set direction of derivation in "dot" belonging spline
//if there are no "next" or "prev" dot set rude derivation
void spdot_t_auto_set_der(spdot_t* dot);

//automatically set direction of derivation in "dot" belonging spline
//if there are no "next" or "prev" do nothing
void spdot_t_auto_set_der_soft(spdot_t* dot);

//return "count" transit spligned points between "cur" and "next" excluding bounding points
p2d* get_transit_spline_points(spdot_t* cur, int count);

//return double[2]* array of 2d points with len = "line_len" consisting spline
double** get_spline(spline_t spline, int* line_len);

#ifdef __cplusplus
}
#endif



#endif
