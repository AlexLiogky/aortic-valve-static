#include <string.h>
#include <assert.h>
#include <stdlib.h>
#include <math.h>
#include <stdio.h>
#include "bound-box.h"

//coord(i, j, k, n) -> (i * n_j * n_k + j * n_k + k) * n_n + n
//coord(n, i, j, k) -> n * n_i * n_j * n_k + i * n_j * n_k + j * n_k + k

void nets_t_fill_borders(nets_t nets, double res[3][2]){
	int i = 0, j = 0, k = 0;
	for (i = 0; i < 3; i++){
		res[i][0] = nets.nets[0].vrtx.nodes[0][0].coord.coord[i];
		res[i][1] = res[i][0];
	}

	int nets_c = nets.count;
	for (i = 0; i < nets_c; i++){
		net_t net = nets.nets[i];
		if (net_is_static(net)) continue;
		int n_cnt = net.vrtx.count;
		for (j = 0; j < n_cnt; j++){
			node_t node = net.vrtx.nodes[j][0];
			for (k = 0; k < 3; k++){
				double coord = node.coord.coord[k];
				if (res[k][0] > coord) res[k][0] = coord;
				if (res[k][1] < coord) res[k][1] = coord;
			}
		}
	}
}

static void box_t_fill_borders(box_t box, double brds[3][2]){
	int i;
	for (i = 0; i < 3; i++){
		brds[i][0] = box.start.coord[i];
		brds[i][1] = box.start.coord[i] + box.step * box.borders[i];
	}
}

static void get_borders(double brds[3][2], box_t box, int res[3]){
	double step = box.step;
	point_t start = box.start;
	int i = 0;
	for (i = 0; i < 3; i++){
		//res[i][0] = (int) floor((brds[i][0] - start.coord[i]) / step);
		res[i] = (int)(ceil((brds[i][1] - start.coord[i]) / step));
		if (res[i] == 0) res[i] = 1;
	}
}

void extend_borders(double coef, double brds[3][2]){
	int i;
	for (i = 0; i < 3; i++){
		brds[i][0] = (brds[i][1] + brds[i][0]) / 2 - coef * (brds[i][1] - brds[i][0]) / 2;
		brds[i][1] = (brds[i][1] + brds[i][0]) / 2 + coef * (brds[i][1] - brds[i][0]) / 2;
	}
}

int nets_t_dyn_cnt(nets_t nets){
	int res = 0;
	for (unsigned int i = 0; i < nets.count; i++)
		res += !net_is_static(nets.nets[i]);
	return res;
}

crd_t box_t_get_crd(box_t box, point_t coord, int net_id){
	crd_t crd;
	point_t shift = point_t_dif(coord, box.start);
	int i;
	for (i = 0; i < 3 ; i++) crd.crd[i] = (int) (floor(shift.coord[i] / box.step));
	crd.crd[3] = net_id;

	//printf("(%d %d %d %d) ", crd.crd[0], crd.crd[1], crd.crd[2], crd.crd[3]);

	return crd;
}

void box_t_get_pnt_ids(box_t box, point_t pnt, int ids[3]){
	point_t shift = point_t_dif(pnt, box.start);
	for (int i = 0; i < 3 ; i++)
		ids[i] = (int) (floor(shift.coord[i] / box.step));
}

void box_t_print_pnt_ids(box_t box, point_t pnt){
	int ids[3];
	box_t_get_pnt_ids(box, pnt, ids);
	printf("point ids = (%d %d %d)\n", ids[0], ids[1], ids[2]);
}

point_t box_t_get_end(box_t box){
	point_t end;
	for (int i = 0; i < 3; i++)
		end.coord[i] = box.start.coord[i] + (box.borders[i] - 1) * box.step;
	return end;
}

int point_belong_box(point_t pnt, box_t box){
	return GEQ(pnt, box.start) && LEQ(pnt, box_t_get_end(box));
}

int crd_is_acceptable(box_t box, crd_t crd){
	int res = 1, i;
	for (i = 0; i < 3; i++)
		res *= ((crd.crd[i] >= 0) && (crd.crd[i] < box.borders[i]));
	res *= ((crd.crd[3] >= 0) && (crd.crd[3] < box.nets_cnt));
	return res;
}

int crd_to_coord(box_t box, crd_t crd){
	int coef[4], res = 0, i;
	coef[0] = 1;
	coef[1] = box.borders[0];
	coef[2] = box.borders[1];
	coef[3] = box.borders[2];
	for (i = 0; i < 4; i++)
		res = res * coef[i] + crd.crd[(i + 3) % 4];
	return res;
}

void init_static(nets_t nets, box_t* box){
	int i, nets_c = nets.count;
	for (i = 0; i < nets_c; i++){
		if (!net_is_static(nets.nets[i])) continue;
		int e_cnt = nets.nets[i].elems.count, j;
		for (j = 0; j < e_cnt; j++){
			point_t n_coord = nets.nets[i].elems.elems[j][0].cntr_mass;
			crd_t crd = box_t_get_crd(box[0], n_coord, i);
			if (crd_is_acceptable(box[0], crd)){
				int b_coord = crd_to_coord(box[0], crd);
				add_elem_to_data(box[0].data, nets.nets[i].elems.elems[j][0].id, b_coord);
			}
		}
	}
}

box_t box_t_construct(nets_t nets, double step){
	double brds[3][2];
	nets_t_fill_borders(nets, brds);
	double ext_coef = 1.1;
	extend_borders(ext_coef, brds);

	box_t box;
	point_t start = point_t_get_point(brds[0][0], brds[1][0], brds[2][0]);
	box.start = start;
	box.step = step;
	box.nets_cnt = nets.count;
	box.dyn_nets_cnt = nets_t_dyn_cnt(nets);
	get_borders(brds, box, box.borders);
	int volume = 1, i;
	for (i = 0; i < 3; i++)
		volume *= box.borders[i];
	box.dyn_data_size = volume * box.dyn_nets_cnt;
	int data_size = volume * nets.count;
	box.data = data_t_construct(data_size, DEFOLT_DATA_ELEM_SIZE);
	init_static(nets, &box);

	return box;
}

void box_t_destruct(box_t* box){
	data_t_destruct((*box).data);
	(*box).step = -1;
	(*box).nets_cnt = -1;
}

/*int box_t_get_coord(box_t box, point_t coord, int net_id){
	point_t shift = point_t_dif(coord, box.start);
	int i, res = (int) floor(shift.coord[0] / box.step);
	for (i = 1; i < 3 ; i++)
		res = (int) floor(shift.coord[i] / box.step) + box.borders[i] * res;
	res = box.nets_cnt * res + net_id;

	return res;
}*/

void recreate_box(nets_t nets, box_t* box, crd_t crd);
void crd_t_dump(crd_t crd);
void box_t_dump(box_t box);

void nets_t_update_box(nets_t nets, box_t* box){
	assert(box != NULL);
	memset(box[0].data.busy_len, 0, box[0].dyn_data_size * sizeof(int));

	int i, nets_c = nets.count;
	for (i = 0; i < nets_c; i++){
		if (net_is_static(nets.nets[i])) continue;
		int e_cnt = nets.nets[i].elems.count, j;
		for (j = 0; j < e_cnt; j++){
			point_t n_coord = nets.nets[i].elems.elems[j][0].cntr_mass;
			crd_t crd = box_t_get_crd(box[0], n_coord, i);
			if (crd_is_acceptable(box[0], crd)){
				int b_coord = crd_to_coord(box[0], crd);
				add_elem_to_data(box[0].data, nets.nets[i].elems.elems[j][0].id, b_coord);
			}
			else{
				recreate_box(nets, box, crd);
				nets_t_update_box(nets, box);
				return;
			}
		}
	}
}

static void triangle_fill_borders(point_t triangle[3], double res[3][2]){
	for (int i = 0; i < 3; i++){
		res[i][0] = triangle[0].coord[i];
		res[i][1] = res[i][0];
	}

	for (int j = 1; j < 3; j++)
		for (int i = 0; i < 3; i++){
			res[i][0] = (res[i][0] < triangle[j].coord[i]) ? res[i][0] : triangle[j].coord[i];
			res[i][1] = (res[i][1] < triangle[j].coord[i]) ? triangle[j].coord[i] : res[i][1];
		}
	return;
}

point_t get_unit_node(point_t unit, point_t step, int x, int y, int z){
	point_t res = unit;
	res.coord[0] += x * step.coord[0];
	res.coord[1] += y * step.coord[1];
	res.coord[2] += z * step.coord[2];
	return res;
}

int _Show = 0;

#define NODE(X, Y, Z) get_unit_node(unit, step, X, Y, Z)
#define SHIFT(X, Y, Z) get_unit_node(get_zero_point(), step, X, Y, Z)
#define INT(NX, NY, NZ, SX, SY, SZ) line_fragment_to_trangle_intersection(NODE(NX, NY, NZ), SHIFT(SX, SY, SZ), triangle, NULL)
int box_unit_intersect_triangle(point_t unit, point_t step, point_t triangle[5]){// triangle[5] = {triangle[3], normal, cntr_mass}
	for (int i = 0; i < 5; i++){
		int flag = 1;
		for (int j = 0; j < 3; j++){
			double dif = triangle[i].coord[j] - unit.coord[j];
			flag *= (dif >= 0) && (dif < step.coord[j]);
		}
		if (flag) return 1;
		if (i == 2) i++;
	}

	//point_t shift = {};

	if (INT(0, 0, 0, 0, 0, 1) || INT(0, 0, 0, 0, 1, 0) || INT(0, 0, 0, 1, 0, 0))
		return 1;
	if ( INT(0, 0, 1, 0, 1, 0) || INT(0, 0, 1, 1, 0, 0)) return 1;
	if (INT(0, 1, 0, 0, 0, 1) ||  INT(0, 1, 0, 1, 0, 0)) return 1;
	if (INT(1, 0, 0, 0, 0, 1) || INT(1, 0, 0, 0, 1, 0)) return 1;
	if (INT(0, 1, 1, 1, 0, 0) || INT(1, 1, 0, 0, 0, 1) || INT(1, 0, 1, 0, 1, 0))
		return 1;

	return 0;
}
#undef NODE
#undef INT

point_t box_t_get_unit_from_crd(box_t box, crd_t crd){
	point_t unit = box.start;
	for (int i = 0; i < 3; i++)
		unit.coord[i] += crd.crd[i] * box.step;
	return unit;
}

#define FOR_DIM(X) for (i[X] = crd_i.crd[X]; i[X] <= crd_f.crd[X]; i[X]++)
int set_exact_elem_to_box(box_t* box, int bar, elem_t elem){ //не проверяет, принадлежит ли треугольник всему боксу
	point_t triangle[5];
	elem_to_triangle(elem, triangle);
	triangle[4] = elem.cntr_mass;
	triangle[3] = elem.or_area;
	double res[3][2];
	triangle_fill_borders(triangle, res);
	point_t start = {{res[0][0], res[1][0], res[2][0]}};
	point_t end = {{res[0][1], res[1][1], res[2][1]}};
	crd_t crd_i = box_t_get_crd(box[0], start, bar);
	crd_t crd_f = box_t_get_crd(box[0], end, bar);
	crd_t curcrd = {};
	curcrd.crd[3] = bar;
	point_t step3d = {{box[0].step, box[0].step, box[0].step}};
	int i[3];
	FOR_DIM(0)
	FOR_DIM(1)
	FOR_DIM(2){
		for (int j = 0; j < 3; j++)
			curcrd.crd[j] = i[j];
		point_t unit = box_t_get_unit_from_crd(*box, curcrd);
		if (box_unit_intersect_triangle(unit, step3d, triangle)){
			int b_coord = crd_to_coord(box[0], curcrd);
			add_elem_to_data(box[0].data, elem.id, b_coord);
		}
	}

	return 1;
}

void box_t_set_elem_exact(nets_t nets, box_t* box){
	assert(box != NULL);
	memset(box[0].data.busy_len, 0, box[0].data.data_cnt * sizeof(int));

	int i, nets_c = nets.count;
	for (i = 0; i < nets_c; i++){
		int e_cnt = nets.nets[i].elems.count, j;
		for (j = 0; j < e_cnt; j++){
			set_exact_elem_to_box(box, i, nets.nets[i].elems.elems[j][0]);
		}
	}
}

void print_brds(double brds[3][2]){
	printf("((%lg %lg) (%lg %lg) (%lg %lg))", brds[0][0], brds[0][1], brds[1][0], brds[1][1], brds[2][0], brds[2][1]);
}

void get_new_borders(nets_t nets, box_t box, double brds[3][2], crd_t crd){
	double newbrds[3][2];
	nets_t_fill_borders(nets, newbrds);
	double oldbrds[3][2];
	box_t_fill_borders(box, oldbrds);
	int i;
	double ext_coef = 1.05;
	for (i = 0; i < 3; i++){
		double newb[2];
		newb[1] = (oldbrds[i][1] + oldbrds[i][0]) / 2 - ext_coef * (oldbrds[i][1] - newbrds[i][0]) / 2;
		newb[0] = (oldbrds[i][1] + oldbrds[i][0]) / 2 - ext_coef * (oldbrds[i][1] - oldbrds[i][0]) / 2;
		brds[i][0] = (oldbrds[i][0] < newbrds[i][0]) ? ((crd.crd[i] < 0) ? newb[0] : oldbrds[i][0]) : newb[1];
		newb[1] = (oldbrds[i][1] + oldbrds[i][0]) / 2 - ext_coef * (oldbrds[i][0] - newbrds[i][1]) / 2;
		newb[0] = (oldbrds[i][1] + oldbrds[i][0]) / 2 - ext_coef * (oldbrds[i][0] - oldbrds[i][1]) / 2;
		brds[i][1] = (oldbrds[i][1] > newbrds[i][1]) ? ((crd.crd[i] >= box.borders[i]) ? newb[1] : oldbrds[i][1]) : newb[1];
	}
}

void recreate_box(nets_t nets, box_t* box, crd_t crd){
	double brds[3][2];
	get_new_borders(nets, box[0], brds, crd);
	point_t start = point_t_get_point(brds[0][0], brds[1][0], brds[2][0]);
	double step = box[0].step;
	int dyn_nets_cnt = box[0].dyn_nets_cnt;
	box_t_destruct(box);

	box[0].start = start;
	box[0].step = step;
	get_borders(brds, box[0], box[0].borders);

	int volume = 1, i;
	for (i = 0; i < 3; i++)
		volume *= box[0].borders[i];
	int data_size = volume * nets.count;
	box[0].dyn_nets_cnt = dyn_nets_cnt;
	box[0].dyn_data_size = volume * dyn_nets_cnt;
	box[0].nets_cnt = nets.count;
	box[0].data = data_t_construct(data_size, DEFOLT_DATA_ELEM_SIZE);
	init_static(nets, box);
}

void crd_t_dump(crd_t crd){
	printf("crd (%d %d %d %d)\n", crd.crd[0], crd.crd[1], crd.crd[2], crd.crd[3]);
}

void box_t_dump(box_t box){
	printf("box:{\n");
	printf("start: "); point_t_dump(box.start);
	printf("step: %lg\n", box.step);
	printf("nets_cnt = %d\n", box.nets_cnt);
	printf("dyn_nets_cnt = %d\n", box.dyn_nets_cnt);
	printf("dyn_data_size = %d\n", box.dyn_data_size);
	printf("borders: ");
	printf("x: (0 %d) "  , box.borders[0]);
	printf("y: (0 %d) "  , box.borders[1]);
	printf("z: (0 %d) \n", box.borders[2]);
	data_t_dump(box.data);
	printf("}\n");
}


#define CHECK_ACCEPTABLE(X, ID)								\
curcrd.crd[X] = pnt_crd.crd[X] + ID;						\
if (!(curcrd.crd[X] >= 0 && curcrd.crd[X] < box.borders[X]))\
	continue

void box_t_local_point_to_net_projection(box_t box, point_t point, net_t net, int net_id, double* sqr_distance){
	crd_t pnt_crd = box_t_get_crd(box, point, net_id), curcrd;
	curcrd.crd[3] = net_id;
	double dist = 1e+20;
	for (int i = -1; i <= 1; i++){ CHECK_ACCEPTABLE(0, i);
	for (int j = -1; j <= 1; j++){ CHECK_ACCEPTABLE(1, j);
	for (int k = -1; k <= 1; k++){ CHECK_ACCEPTABLE(2, k);
		int coord = crd_to_coord(box, curcrd);
		unsigned int e_cnt = box.data.busy_len[coord];
		double loc_dist = 0;
		if (e_cnt != 0) point_to_elem_projection(point, net.elems.elems[box.data.data[coord][0]], &loc_dist);
		else continue;
		for (unsigned int l = 1; l < e_cnt; l++){
			double cur_dist = 1e+25;
			point_to_elem_projection(point, net.elems.elems[box.data.data[coord][l]], &cur_dist);
			if (cur_dist < loc_dist) loc_dist = cur_dist;
		}
		if (e_cnt != 0 && loc_dist < dist) dist = loc_dist;
	}}}
	if (sqr_distance != NULL) *sqr_distance = dist;
}
#undef CHECK_ACCEPTABLE

void box_t_local_node_to_net_projection(box_t box, node_t* node, net_t net, int net_id, double* sqr_distance){
	box_t_local_point_to_net_projection(box, node[0].coord, net, net_id, sqr_distance);
}

