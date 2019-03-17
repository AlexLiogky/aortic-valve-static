#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <math.h>
#include "precomputation.h"
#include "intersection.h"
#include "sew_leaf.h"

//##############################dumpers#################################

void segline_t_dump(segline_t* line){
	if (!line) {printf("segline {NULL}\n"); return;}
	printf("segline {\n");
	printf("from: "), point_t_dump(line[0].rep);
	printf("len = %d, &line = %p\n", line[0].len, line[0].line);
	if (line[0].updated) printf("status: updated\n");
	else printf("status: outdated\n");
	printf("next: %p\n", line[0].next);
	printf("}\n");
}

void segline_t_full_dump(segline_t* line){
	if (!line) {printf("segline: {NULL}\n"); return;}
	printf("segline_t:{ \n");
	segline_t* next = NULL; 
	if (line) next = line[0].next;
	printf("head:: ");
	segline_t_dump(line);
	while (next){
		printf("next:: ");
		segline_t_dump(next);
		next = next[0].next;
	}
	printf("}\n");
}

void sew_line_t_plate_params_dump(sew_line_t* sew){
	if (!sew || !sew[0].curvilinear_coords_params) 
		{printf("params: NULL\n"); return;}
	point_t start = ((point_t*) (sew[0].curvilinear_coords_params))[0];
	point_t base[3];
	for (int i = 1; i <= 3; i++)
		base[i-1] = ((point_t*) sew[0].curvilinear_coords_params)[i];
	printf("start point: "); point_t_dump(start);
	for (int i = 1; i <= 3; i++)
		printf("base[%d]: ", i), point_t_dump(base[i-1]);
} 

void sew_line_t_dump(sew_line_t* sew, void (*params_dump)(struct sew_line_t* sew)){
	if (!sew) {printf("sew_line_t {NULL}\n"); return;}
	printf("sew_line_t {\n");
	printf("aorta "); nets_t_dump(sew[0].aorta);
	printf("box "); box_t_dump(sew[0].aorta_box);
	if (params_dump) {
		printf("params:{\n");
		params_dump(sew);
		printf("}\n");
	}
	segline_t_full_dump(sew[0].head);
	spline_t_full_dump(sew[0].sew2d);
	printf("}\n");
}

void sew_line_t_plate_dump(sew_line_t* sew){
	sew_line_t_dump(sew, sew_line_t_plate_params_dump);
}

//###########################dumpers####################################

segline_t* segline_t_construct(unsigned int len, point_t rep, point_t* line, segline_t* next, int update){
	segline_t* sline = (segline_t*)malloc(sizeof(segline_t));
	sline[0].rep = rep;
	sline[0].len = len;
	sline[0].line = line;
	sline[0].next = next;
	sline[0].updated = update;
	return sline;
}

void segline_t_update(segline_t* segline, unsigned int len, point_t* line){
	if (len < 1) return;
	if (segline[0].line) free(segline[0].line);
	segline[0].len = len;
	segline[0].line = line;
	segline[0].updated = 1;
}

void segline_t_set_outdate(segline_t* line){
	if (line) line[0].updated = 0;
}

void add_segline_t(segline_t* newseg, segline_t* prev, segline_t* prevprev){
	assert(prev);
	newseg[0].next = prev[0].next;
	prev[0].next = newseg;
	segline_t_set_outdate(prevprev);
	segline_t_set_outdate(prev);
	segline_t_set_outdate(newseg);
	segline_t_set_outdate(newseg[0].next);
}

line_t segline_t_get_line(segline_t* head){
	int len = 0;
	segline_t* cur = head;
	while (cur){
		len += 1 + cur[0].len;
		cur = cur[0].next;
	}
	line_t line = {len, (point_t*)calloc(len, sizeof(point_t))};
	
	cur = head;
	int off = 0;
	while (cur){
		line.line[off++] = cur[0].rep;
		for (int i = 0; i < cur[0].len; i++)
			line.line[off++] = cur[0].line[i];
		cur = cur[0].next;
	}
	
	return line;
}


p2d convert_3d_to_2d_plate(point_t p, sew_line_t* sew){
	point_t start = ((point_t*) (sew[0].curvilinear_coords_params))[0];
	point_t base[2];
	for (int i = 1; i <= 2; i++)
		base[i-1] = ((point_t*) sew[0].curvilinear_coords_params)[i];
	point_t rel = DIF(p, start);
	p2d res = {{DOT(rel, base[0]), DOT(rel, base[1])}};
	return res;
}

point_t convert_2d_to_3d_plate(p2d p, sew_line_t* sew){
	point_t start = ((point_t*) sew[0].curvilinear_coords_params)[0];
	point_t base[2];
	for (int i = 1; i <= 2; i++)
		base[i-1] = ((point_t*) sew[0].curvilinear_coords_params)[i];
	point_t res = SUM(SUM(SCAL(p.pnt[0], base[0]), SCAL(p.pnt[1], base[1])), start);
	return res;
}

int place_point_to_aorta_plate(point_t p, sew_line_t* sew, point_t* aorta_p){
	point_t line[2];
	line[0] = p;
	line[1] = ((point_t*) sew[0].curvilinear_coords_params)[3];
	line[1] = SUM(line[0], line[1]);
	point_t swap = line[0];
	line[0] = line[1];
	line[1] = swap;
		
	if (!line_to_boxed_nets_intersection(line, sew[0].aorta, 0, sew[0].aorta_box, aorta_p)){
		printf("Intersection doesn't detected\n");
		return 0;
	}
	/*static int a = 0;
	static double max_sqr_cntr_mass_dist = 0;
	if (a == 0) max_sqr_cntr_mass_dist = get_max_sqr_cntr_mass_dist(sew[0].aorta.nets[0]);
	a++;
	point_t swap = line[0];
	line[0] = line[1];
	line[1] = swap;
	*aorta_p = find_line_to_net_intersection(sew[0].aorta.nets[0], line, max_sqr_cntr_mass_dist, NULL);
	{printf("Intersection: ");
	point_t_dump(line[0]);
	point_t_dump(line[1]);
	point_t_dump(*aorta_p);}*/
	
	return 1;
}

sew_line_t construct_plate_approx(point_t start, point_t base[3], point_t ends[2], point_t ends_der[2]){	//base should be orthonormal
	sew_line_t sew;
	sew.curvilinear_coords_params = malloc(4 * sizeof(point_t));
	((point_t*) sew.curvilinear_coords_params)[0] = start;
	for (int i = 1; i <= 3; i++)
		((point_t*) sew.curvilinear_coords_params)[i] = base[i-1];
	sew.convert_3d_to_2d = &convert_3d_to_2d_plate;
	sew.convert_2d_to_3d = &convert_2d_to_3d_plate;
	sew.place_point_to_aorta = &place_point_to_aorta_plate;
	
	p2d pnt = sew.convert_3d_to_2d(ends[0], &sew);
	p2d der = sew.convert_3d_to_2d(SUM(ends_der[0], start), &sew);
	spdot_t* first = spdot_t_construct(pnt.pnt, der.pnt);
	sew.sew2d = spline_construct(first);
	pnt = sew.convert_3d_to_2d(ends[1], &sew);
	der = sew.convert_3d_to_2d(SUM(ends_der[1], start), &sew);
	spdot_t* end = spdot_t_construct(pnt.pnt, der.pnt);
	spline_t_add_dot(&sew.sew2d, end, first);
	segline_t* segend = segline_t_construct(0, ends[1], NULL, NULL, 1);
	sew.head = segline_t_construct(0, ends[0], NULL, segend, 0);
	
	return sew;
}


void sewline_t_update_seglines(sew_line_t sew){
	//printf("Hey2\n");
	spdot_t* spcur = sew.sew2d.head;
	while (spcur){
		spdot_t_auto_set_der_soft(spcur);
		spcur = spcur[0].next;
	}
	//printf("Hey3\n");
	segline_t* segcur = sew.head;
	spcur = sew.sew2d.head;
	const int segcnt = 50;
	while (segcur){
		//printf("Hey4\n");
		if (!segcur[0].updated && segcur[0].next){
			//printf("Hey5\n");
			p2d* pnts = get_transit_spline_points(spcur, segcnt);
			/*printf("pnts:: \n");
			for (int i = 0; i < segcnt; i++){
				p2d_dump(pnts[i]);
			} */
			//printf("Hey6, segcnt = %d\n", segcnt);
			point_t* projs = (point_t*)malloc(segcnt * sizeof(point_t));
			//printf("Hey6.5 pnts = %p, projs = %p\n", pnts, projs);
			for (int i = 0; i < segcnt; i++){
				projs[i] = sew.convert_2d_to_3d(pnts[i], &sew);
				if(!sew.place_point_to_aorta(projs[i], &sew, projs + i))
					printf("Point is not projectable\n");
			}
			free(pnts);
			/*printf("projs:: \n");
			for (int i = 0; i < segcnt; i++){
				point_t_dump(projs[i]);
			} */
			//printf("Hey7\n");
			//printf("Hey8\n");
			/*printf("updated projs:: \n");
			for (int i = 0; i < segcnt; i++){
				point_t_dump(projs[i]);
			} */
			segline_t_update(segcur, segcnt, projs);
			//printf("Hey9\n");
		}
		if (!segcur[0].next) segline_t_update(segcur, 0, NULL);
		//printf("Hey10\n");
		segcur = segcur[0].next;
		spcur = spcur[0].next;
	}
	//printf("Hey11\n");
}

void sew_line_t_add_point(point_t p, sew_line_t* sew){
	if (!sew) return;
	p2d new_p = sew[0].convert_3d_to_2d(p, sew);
	int id = 0;
	
	spdot_t* spprev = get_nearest_prev_dot(sew[0].sew2d, new_p.pnt, &id);
	spdot_t* dot = spdot_t_construct(new_p.pnt, NULL);
	spline_t_add_dot(&sew[0].sew2d, dot, spprev);
	
	segline_t* prevprev = NULL;
	segline_t* prev = NULL;
	segline_t* cur = sew[0].head;
	for (int i = 0; i <= id; i++){
		if (i == id - 1) prevprev = cur;
		if (i == id) prev = cur;
		cur = cur[0].next; 
	}
	
	segline_t* newseg = segline_t_construct(1, p, NULL, prev, 0);
	add_segline_t(newseg, prev, prevprev);
}

line_t sew_line_t_get_line_on_aorta(sew_line_t* sew){
	//printf("Hey1\n");
	sewline_t_update_seglines(sew[0]);
	//printf("Hey-inf\n");
	//sew_line_t_plate_dump(&sews[0]);
	
	return segline_t_get_line(sew[0].head);
}

sew_line_t* init_sew_line_on_aorta(nets_t aorta, point_t* commissur, point_t* bottom, int cnt_leafs){
	double* a = (double*)calloc(cnt_leafs, sizeof(double));
	point_t* base = (point_t*)calloc(cnt_leafs, sizeof(point_t)); 
	point_t* normal = (point_t*)calloc(cnt_leafs, sizeof(point_t));
	point_t* complement = (point_t*)calloc(cnt_leafs, sizeof(point_t));
	for (int i = 0; i < cnt_leafs; i++){
		base[i] = DIF(commissur[(i+1) % cnt_leafs], commissur[i % cnt_leafs]);
		a[i] = LEN(base[i]);
		SCAL_S(1.0 / a[i], &(base[i]));
		normal[i] = NORM(OR_AREA(commissur[i % cnt_leafs], commissur[(i + 1) % cnt_leafs], bottom[i]));
		complement[i] = CROSS(normal[i], base[i]);
	}
	sew_line_t* sews = (sew_line_t*)calloc (cnt_leafs, sizeof(sew_line_t));
	
	if (cnt_leafs >= 3){
		point_t axis = OR_AREA(commissur[0], commissur[1], commissur[2]);
		double R = a[0] * a[1] * a[2] / (4 * LEN(axis));
		NORM_S(&axis);
		for (int i = 0; i < cnt_leafs; i++){
			double d = R, b = R / DOT(normal[i], axis), c = a[i] / 2;
			double x0 = c, y0 = -b * sqrt(1 - (c / d) * (c / d));
			double tgf = - (d / b) * (d / b) * y0 / x0;
			point_t inc = SCAL(tgf, base[i]);
			
			point_t start = commissur[i];
			point_t end_dir[2] = {NORM(ADD(inc, -1, complement[i])), NORM(ADD(inc, 1, complement[i]))};
			//for (int j = 0; j < 2; j++)
				//end_dir[j] = SUM(end_dir[j], start);		//for correct converting in 2d point
			printf("normals: \n"); point_t_dump(end_dir[0]); point_t_dump(end_dir[1]);
			point_t basis[3] = {base[i], complement[i], normal[i]};
			point_t ends[2] = {commissur[i % cnt_leafs], commissur[(i + 1) % cnt_leafs]};
			sews[i] = construct_plate_approx(start, basis, ends, end_dir);
		}
	}
	else {
		for (int i = 0; i < cnt_leafs; i++){
			point_t basis[3] = {base[i], complement[i], normal[i]};
			point_t start = commissur[i];
			point_t end_dir[2] = {SCAL(-1, complement[i]), complement[i]};
			//for (int j = 0; j < 2; j++)
				//end_dir[j] = SUM(end_dir[j], start);	//for correct converting in 2d point
			point_t ends[2] = {commissur[i % cnt_leafs], commissur[(i + 1) % cnt_leafs]};
			sews[i] = construct_plate_approx(start, basis, ends, end_dir);
		}
	}
	
	box_t aorta_box = box_t_construct(aorta, get_Contact_Resolution());
	//nets_t_update_box(aorta, &aorta_box);
	box_t_set_elem_exact(aorta, &aorta_box);
	for (int i = 0; i < cnt_leafs; i++){
		sews[i].aorta = aorta;
		sews[i].aorta_box = aorta_box;
		sew_line_t_add_point(bottom[i], &sews[i]);
	}
	
	return sews;
}

//TODO: деструктор

