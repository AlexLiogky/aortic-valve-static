#include <math.h>
#include <float.h>
#include <stdlib.h>
#include "intersection.h"

#define K moving[k]
#define O slc_i
int _find_line_box_borders(int borders[2][2], point_t l_dir, point_t r0,\
int slc_i, int slc, int* moving, box_t box, double* step){
	int res = 1;
	for (int k = 0; k < 2; k++){				
		double base = (l_dir.coord[K]*(slc*step[O] - r0.coord[O]) + r0.coord[K])/step[K];
		double range = l_dir.coord[K] * step[O] / step[K];
		double start = base - 1 + ((range < 0) ? range : 0);
		double end =  base + ((range > 0) ? range : 0);
		borders[k][0] = (int) ceil(start);
		borders[k][1] = (int) floor(end);
		if (borders[k][0] < 0) borders[k][0] = 0;
		if (borders[k][1] >= box.borders[K]) borders[k][1] = box.borders[K] - 1;
		res *= (borders[k][1] >= borders[k][0]);
		}
	return res;
}
#undef K
#undef O

int _check_slice_box_intersection(point_t line[2], int slc, int slc_i, int* moving, int bar,\
point_t l_dir, point_t r0, nets_t nets, box_t box, double* step, int* flag, point_t* intersect){
	//if (slc == 142) printf("HEY\n");
	int borders[2][2] = {};									
	*flag = _find_line_box_borders(borders, l_dir, r0, slc_i, slc, moving, box, step);
	//if (slc == 80) printf ("borders ((%d %d), (%d %d)) slc = %d\n", borders[0][0], borders[0][1], borders[1][0], borders[1][1], slc);
	for (int l = borders[0][0]; l <= borders[0][1]; l++)
	for (int p = borders[1][0]; p <= borders[1][1]; p++){
		crd_t curcrd;
		curcrd.crd[slc_i] = slc;
		curcrd.crd[moving[0]] = l;
		curcrd.crd[moving[1]] = p;
		curcrd.crd[3] = bar;
		int coord = crd_to_coord(box, curcrd);
		data_t data = box.data;
		unsigned int e_cnt = data.busy_len[coord];
		//printf("e_cnt = %d\n", e_cnt);
		//if (slc == 61 && p == 133) printf ("elem cnt = %d\n", e_cnt);
		//if (slc == 80) printf ("elem cnt = %d\n", e_cnt);
		for (unsigned int q = 0; q < e_cnt; q++){
			int elem_id = data.data[coord][q];
			//if (slc == 80) printf("elem_id = %d\n", elem_id);
			elem_t* elem = nets.nets[bar].elems.elems[elem_id];
			if (line_to_elem_intersection(line, elem, intersect))
				return 1;
		}
	}
	return 0;
}

int _get_line_point(point_t st, point_t dir, int max_i, box_t box, point_t* intersect){
	if (point_belong_box(st, box)){
		if (intersect) *intersect = st;
		return 1;
	}
	point_t start = DIF(box.start, st), end = DIF(box_t_get_end(box), st);
	double t_min = start.coord[max_i] / dir.coord[max_i], t_max = end.coord[max_i] / dir.coord[max_i];
	for (int i = 0; i < 3; i++){
		if (fabs(dir.coord[i]) > DBL_EPSILON){
			t_min = fmax(start.coord[i] / dir.coord[i], t_min);
			t_max = fmin(end.coord[i] / dir.coord[i], t_max);
		}
		else if (start.coord[i] > 0 || end.coord[i] < 0)
			return 0;
		if (t_min > t_max) return 0;
	}
	if (intersect) *intersect = ADD(st, t_min, dir);
	return 1;
}

int line_to_boxed_nets_intersection(point_t line[2], nets_t nets, int bar, box_t box, point_t* intersect){
	int max_i = -1;
	double step[3] = {box.step, box.step, box.step};
	//box_t_dump(box);
	point_t l_dir = signed_norm_inf(DIF(line[1], line[0]), NULL, &max_i);
	//printf("l_direction "); point_t_dump(l_dir);
	int moving[2], off = 0;
	for (int i = 0; i < 3; i++)
		if (i != max_i) moving[off++] = i;
	point_t r0;						
	if (!_get_line_point(line[0], l_dir, max_i, box, &r0))
		return 0;
	//printf("start = "); point_t_dump(r0);
	crd_t crd_r0 = box_t_get_crd(box, r0, bar);
	//printf("start crds: "); box_t_print_pnt_ids(box, r0);
	r0 = DIF(r0, box.start);
	//printf("relative start: "); point_t_dump(r0);
	//printf("max_i = %d\n", max_i);
	
	int flag_p = 1, flag_n = 1;
	for (int i_p = crd_r0.crd[max_i] + 1, i_n = crd_r0.crd[max_i];\
					i_p < box.borders[max_i] || i_n >= 0; i_p++, i_n--){
		if (i_p < box.borders[max_i] && flag_p && \
		_check_slice_box_intersection(line, i_p, max_i, moving, bar,\
		l_dir, r0, nets, box, step, &flag_p, intersect))
				return 1;
		
		if (i_n >= 0 && flag_n && \
		_check_slice_box_intersection(line, i_n, max_i, moving, bar,\
		l_dir, r0, nets, box, step, &flag_n, intersect))
				return 1;
		//printf("i_p = %d, i_n = %d\n", i_p, i_n);
	}
	return 0;
}


/*
 * 			int borders[2][2];									
			flag_p = _find_line_box_borders(borders, l_dir, r0, max_i, moving, box);
			for (l = borders[0][0]; l <= borders[0][1]; l++)
			for (p = borders[1][0]; p <= borders[1][1]; p++){
				crd_t curcrd;
				curcrd.crd[max_i] = i_p;
				curcrd.crd[moving[0]] = l;
				curcrd.crd[moving[1]] = p;
				curcrd.crd[3] = bar;
				int coord = crd_to_coord(box, curcrd);
				data_t data = box.data;
				unsigned int e_cnt = data.busy_len[coord];
				for (int q = 0; q < e_cnt; q++){
					int elem_id = data.data[coord][q];
					elem_t* elem = nets.nets[bar].elems.elems[elem_id];
					if (line_to_elem_intersection(line, elem, intersect)) 
						return 1;
				}
			}
*/

