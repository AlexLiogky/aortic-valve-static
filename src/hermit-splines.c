#include <assert.h>
#include <stdlib.h>
#include <stdio.h>
#include "point.h"
#include "float.h"
#include "hermit-splines.h"


#define HDIM 2
#define INIT_MEM 256

void p2d_dump(p2d obj){
	printf("(%lg; %lg)\n", obj.pnt[0], obj.pnt[1]);
}

void spdot_t_dump(spdot_t* obj){
	printf("spdot_t: %p {\n", obj);
	if (obj){
		printf("pnt ");
		p2d_dump(obj[0].pnt);
		printf("der ");
		p2d_dump(obj[0].der);
		printf("next: %p\n", obj[0].next);
		printf("prev: %p\n", obj[0].prev);
	}
	printf("}\n");
}

void spline_t_full_dump(spline_t sp){
	printf("spline_t {\n");
	printf("len = %d, head = %p\n", sp.len, sp.head);
	spdot_t* dot = sp.head;
	while (dot){
		spdot_t_dump(dot);
		dot = dot[0].next;
	}
	printf("}\n");
}

spline_t spline_construct(spdot_t* head){
	spline_t spline = {1, head};
	return spline;
}

void spdot_destruct(spdot_t* dot){
	free(dot);
}

p2d get_p2d(double pnt[HDIM]){
	p2d p;
	for (int i = 0; i < HDIM; i++)
		p.pnt[i] = pnt[i];
	return p;
}

spdot_t* spdot_t_construct(double pnt[HDIM], double der[HDIM]){
	assert(pnt);
	spdot_t* dot = (spdot_t*)calloc(sizeof(spdot_t), 1);
	(*dot).pnt = get_p2d(pnt);
	if (der) (*dot).der = get_p2d(der);
	(*dot).prev = NULL;
	(*dot).next = NULL;
	return dot;
}

void spdot_t_set_der(spdot_t* dot, double der[HDIM]){
	dot[0].der = get_p2d(der);
}

void spline_destruct(spline_t spline){
	spdot_t* cur = spline.head;
	do{
		spdot_t* next = cur[0].next;
		spdot_destruct(cur);
		cur = next;	
	} while (cur);
	spline.len = -1;
}

void spline_t_add_dot(spline_t* spline, spdot_t* insert, spdot_t* prev){
	assert(insert && prev);
	spline[0].len++;
	insert[0].prev = prev;
	insert[0].next = prev[0].next;
	prev[0].next = insert;
	if (insert[0].next)
		insert[0].next[0].prev = insert;
}

point_t extend_2d_to_3d_point(p2d p){
	return point_t_get_point(p.pnt[0], p.pnt[1], 0);
}

p2d cut_3d_to_2d_point(point_t p){
	return get_p2d(p.coord);
}

#define CUT(X) cut_3d_to_2d_point(X)
#define EXT(X) extend_2d_to_3d_point(X)
double sqr_dist_to_2D_line_fragment(p2d frag[2], p2d p){
	point_t pnt = EXT(p);
	point_t line_segment[2] = {EXT(frag[0]), EXT(frag[1])};
	point_t proj = point_to_line_segment_projection(pnt, line_segment);
	return SQR_LEN(DIF(pnt, proj));
}


spdot_t* get_nearest_prev_dot(spline_t spline, double pnt[HDIM], int* id){
	if (spline.len <= 2) {
		if (id) *id = 0;
		return spline.head;
	}
	
	spdot_t* res = spline.head;
	spdot_t* cur = spline.head;
	spdot_t* next = cur[0].next;
	p2d p = get_p2d(pnt);
	p2d frag[2] = {cur[0].pnt, next[0].pnt};
	double dist = sqr_dist_to_2D_line_fragment(frag, p);
	int spdot_id = 0, k = 0;
	do{
		frag[0] = cur[0].pnt;
		frag[1] = next[0].pnt;
		double curdist = sqr_dist_to_2D_line_fragment(frag, p);
		if (curdist - dist < DBL_EPSILON){
			dist = curdist;
			res = cur;
			spdot_id = k;
		}
		cur = next;
		next = cur[0].next;	
		k++;
	} while (next);
	
	if (id) *id = spdot_id;
	
	return res;
}

void spdot_t_auto_set_der(spdot_t* dot){
	p2d cur = dot[0].pnt;
	if (dot[0].next == NULL){
		if (dot[0].prev == NULL){
			double der[HDIM] = {0, 1};
			dot[0].der = get_p2d(der);
			return;
		}
		p2d prev = dot[0].prev[0].pnt;
		p2d der = CUT(NORM(DIF(EXT(cur), EXT(prev))));
		dot[0].der = der;
		return;
	}
	else if (dot[0].prev == NULL){
		p2d next = dot[0].next[0].pnt;
		p2d der = CUT(NORM(DIF(EXT(next), EXT(cur))));
		dot[0].der = der;
		return;
	}
	point_t next = DIF(EXT(dot[0].next[0].pnt), EXT(cur));
	point_t prev = DIF(EXT(cur), EXT(dot[0].prev[0].pnt));
	double next_len = LEN(next), prev_len = LEN(prev);
	double len = next_len + prev_len;
	SCAL_S(next_len / len, &next);
	SCAL_S(prev_len / len, &prev);
	p2d der = CUT(SUM(next, prev));
	dot[0].der = der;
	return;
}

void spdot_t_auto_set_der_soft(spdot_t* dot){
	if ((dot[0].prev == NULL || dot[0].next == NULL) && dot[0].der.pnt)
		return;
	spdot_t_auto_set_der(dot);
}

p2d* get_transit_spline_points(spdot_t* cur, int count){
	p2d pnt1 = cur[0].pnt, der1 = cur[0].der;
	spdot_t* next = cur[0].next;
	if (!next) return NULL;
	p2d* res = (p2d*)malloc(count * sizeof(p2d));
	p2d pnt2 = next[0].pnt, der2 = next[0].der;
	double step = 1.0 / (count + 1);
	for (int i = 1; i <= count; i++){
		double w = step * i;
		double g0 = 1 - 3*w*w + 2*w*w*w;
		double g1 = w - 2*w*w + w*w*w;
		double h0 = 3*w*w - 2*w*w*w;
		double h1 = -w*w + w*w*w;
		for (int j = 0; j < HDIM; j++){
			res[i-1].pnt[j] = g0 * pnt1.pnt[j] + g1 * der1.pnt[j] +\
			h0 * pnt2.pnt[j] + h1 * der2.pnt[j];
		}
	}
	
	return res;
}

void _set_transit_pnts(spdot_t* cur, int freq, double* pnts[HDIM]){
	p2d pnt1 = cur[0].pnt, der1 = cur[0].der;
	spdot_t* next = cur[0].next;
	if (!next) return;
	p2d pnt2 = next[0].pnt, der2 = next[0].der;
	double step = 1.0 / (freq + 1);
	for (int i = 1; i <= freq; i++){
		double w = step * i;
		double g0 = 1 - 3*w*w + 2*w*w*w;
		double g1 = w - 2*w*w + w*w*w;
		double h0 = 3*w*w - 2*w*w*w;
		double h1 = -w*w + w*w*w;
		for (int j = 0; j < HDIM; j++){
			pnts[i-1][j] = g0 * pnt1.pnt[j] + g1 * der1.pnt[j] +\
			h0 * pnt2.pnt[j] + h1 * der2.pnt[j];
		}
	}
}

double** get_spline(spline_t spline, int* line_len){
	assert(line_len);
	if (spline.len < 2){
		*line_len = 1;
		double** pnts = (double**)malloc(sizeof(double[HDIM]));
		pnts[0] = spline.head[0].pnt.pnt;
		return pnts;
	}
	int freq = *line_len / (spline.len - 1);
	*line_len = freq * (spline.len - 1) + spline.len;
	double** pnts = (double**)malloc((*line_len)*sizeof(double*));
	for (int i = 0; i < *line_len; i++)
		pnts[i] = (double*)malloc(sizeof(double[HDIM]));
	
	spdot_t* cur = spline.head;
	int off = 0;
	for (int j = 0; j < HDIM; j++) pnts[off][j] = cur[0].pnt.pnt[j];
	off++;
	for (int i = 0; i < spline.len - 1; i++){
		_set_transit_pnts(cur, freq, &pnts[off]);
		cur = cur[0].next;
		off += freq;
		for (int j = 0; j < HDIM; j++) pnts[off][j] = cur[0].pnt.pnt[j];
		off++;
	}
	
	return pnts;
}

#undef CUT
#undef EXT
#undef HDIM
#undef INIT_MEM

