#ifndef _SEWLEAF_H
#define _SEWLEAF_H

#include "hermit-splines.h"
#include "nets.h"
#include "bound-box.h"

#ifdef __cplusplus
extern "C" {
#endif

typedef struct segline_t{
	int len;
	point_t rep;
	point_t* line;
	struct segline_t* next;
	int updated;
}segline_t;

typedef struct sew_line_t{
	nets_t aorta;
	box_t aorta_box;
	spline_t sew2d;
	void* curvilinear_coords_params;
	p2d (*convert_3d_to_2d)(point_t p, struct sew_line_t* sew);
	point_t (*convert_2d_to_3d)(p2d p, struct sew_line_t* sew);
	int (*place_point_to_aorta)(point_t p, struct sew_line_t* sew, point_t* aorta_p);
	segline_t* head;
}sew_line_t;

typedef struct sews_t{
	int count;
	sew_line_t* sew;
}sews_t;

//debug print of sew_line
//params_dump - function for printing some sew information
void sew_line_t_dump(sew_line_t* sew, void (*params_dump)(struct sew_line_t* sew));
void sew_line_t_plate_params_dump(sew_line_t* sew); //debug print of sew plate parameters
void sew_line_t_plate_dump(sew_line_t* sew); //debug print of sew_line with plate information
void segline_t_dump(segline_t* line); //debug print of segline_t
void segline_t_full_dump(segline_t* line); //debug detailed print of segline_t

//creates dynamicaly segline
//len - length of line
//rep - reference point of segline
//line - array of points forming line
//next - pointer to next segment of sewed line
//update - status flag
segline_t* segline_t_construct(unsigned int len, point_t rep, point_t* line, segline_t* next, int update);

//set "len", "line" and change update to true
void segline_t_update(segline_t* segline, unsigned int len, point_t* line);

//check "update" status
void segline_t_set_outdate(segline_t* line);

//add newseg after prev and set corresponding status where neccessary
void add_segline_t(segline_t* newseg, segline_t* prev, segline_t* prevprev);

//return line stored in segline's list starting from "head"
line_t segline_t_get_line(segline_t* head);

//add point "p" on aorta to sew line "sew" to fit result
void sew_line_t_add_point(point_t p, sew_line_t* sew);

//return line of sewing on aorta of corresponding leaflet
line_t sew_line_t_get_line_on_aorta(sew_line_t* sew);

//return initialized array sew_line_t[cnt_leafs]
//aorta - only nets containing aorta
//commissur - points of commissure on aorta
//bottom - lower sewing points placed on commissur[i] - commissur[i+1] sewed line
//cnt_leafs - count of stitched leaflets
sew_line_t* init_sew_line_on_aorta(nets_t aorta, point_t* commissur, point_t* bottom, int cnt_leafs);

//TODO: деструктор

#ifdef __cplusplus
}
#endif


#endif
