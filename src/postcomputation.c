#include <stdlib.h>
#include <math.h>
#include <stdio.h>
#include <float.h>
#include "precomputation.h"
#include "postcomputation.h"
#include "bound-box.h"

double node_t_area(net_t net, node_t* node){
	unsigned int e_cnt = (*node).cnt_elems, i;
	double area = 0;
	for (i = 0; i < e_cnt; i++){
		area = point_t_length((*(net.elems.elems[(*node).elems_id[i]])).or_area);
	}
	return area / 3;
}

int elem_t_check_coapt(elem_t* elem){
	unsigned int check = 1, i;
	for (i = 0; i < 3; i++)
		check *= ((*((*elem).vrts[i])).coaptative + 1);
	return (check > 0);
}

void print_statistic(nets_t nets){
	double mid_deform = 0;
	spring_t* some_spr = nets.nets[0].springs.springs[0];
	double some_deform = spring_t_get_deformation(some_spr);
	double min_deform = some_deform, max_deform = some_deform;
	int cnt = 0;
	for (unsigned int i = 0; i < nets.count; i++)
		for (unsigned int j = 0; j < nets.nets[i].springs.count; j++){
			spring_t* spr = nets.nets[i].springs.springs[j];
			if (is_fix(spr[0].ends[0][0].state) && is_fix(spr[0].ends[1][0].state)) continue;
			double deform = spring_t_get_deformation(spr);
			min_deform = (min_deform > deform)? deform : min_deform;
			max_deform = (max_deform < deform)? deform : max_deform;
			//if (deform < 0) printf("i = %d, d = %lg\n", i,  deform);
			mid_deform += fabs(deform);
			cnt++;
		}
	mid_deform /= cnt;
	printf("min_deformation = %lg\n", min_deform);
	printf("max_deformation = %lg\n", max_deform);
	printf("mid_deformation = %lg\n", mid_deform);
}
//#############################################################################
void node_t_set_coaptative(node_t* node, unsigned int bar, unsigned int nets_cnt){
	if (node[0].coaptative < 0) node[0].coaptative = 0;
	node[0].coaptative = node[0].coaptative * nets_cnt + bar;
}

int node_t_check_coaptative_with(node_t* node, unsigned int bar, unsigned int nets_cnt){
	int coapt_code = node[0].coaptative;
	int coapt = 0;
	if (coapt_code < 0) return coapt;
	int coapt_num = 0;
	do{
		coapt += (coapt_code % nets_cnt == bar);
		coapt_num++;
		coapt_code /= nets_cnt;
	} while (coapt_code > 0);

	return coapt * coapt_num;
}

int elem_t_check_coapt_with(elem_t* elem, int bar, int nets_cnt){
	int check = 1;
	for (int i = 0; i < 3; i++)
		check *=  node_t_check_coaptative_with((*elem).vrts[i], bar, nets_cnt);
	return (check > 0);
}

int node_t_get_coaptative_power(node_t* node, unsigned int nets_cnt){
	int coapt_code = node[0].coaptative;
	int coapt_pow = 0;
	if (coapt_code < 0) return coapt_pow;
	do{
		coapt_pow++;
		coapt_code /= nets_cnt;
	} while (coapt_code > 0);

	return coapt_pow;
}

void nets_t_set_coaptative(nets_t nets){				//нужно изменить node_to_net_projection с учётом box_t
	box_t box = box_t_construct(nets, get_Contact_Resolution());
	nets_t_update_box(nets, &box);

	unsigned int nets_cnt = nets.count;
	if (nets_cnt < 2) return;
	for (unsigned int i = 0; i < nets_cnt; i++){
		net_t origin = nets.nets[i];
		unsigned int node_cnt = origin.vrtx.count, k;
		for (k = 0; k < node_cnt; k++){
			node_t* node = origin.vrtx.nodes[k];
			node[0].coaptative = -1;
		}
	}

	for (unsigned int i = 0; i < nets_cnt; i++)
		for (int j = (int)nets_cnt - 1; j >= 0; j--){
			if (j == (int)i) continue;
			net_t origin = nets.nets[i];
			net_t barrier = nets.nets[j];
			unsigned int node_cnt = origin.vrtx.count;
			for (unsigned int k = 0; k < node_cnt; k++){
				node_t* node = origin.vrtx.nodes[k];
				double dist = 2e+20;		//very big number
				box_t_local_node_to_net_projection(box, node, barrier, j, &dist);
				if (dist < Contact_Resolution) node_t_set_coaptative(node, j, nets_cnt);
			}
		}

	box_t_destruct(&box);
	return;
}

vrtx_t get_coapt_field(nets_t nets, int cur, int bar){
	net_t net = nets.nets[cur];
	node_t** nodes = net.vrtx.nodes;
	unsigned int n_cnt = net.vrtx.count, i, count = 0;
	for (i = 0; i < n_cnt; i++)
		if (node_t_check_coaptative_with(nodes[i], bar, nets.count)) count++;
	node_t** field_n = (node_t**)calloc(count, sizeof(node_t*));
	count = 0;
	for (i = 0; i < n_cnt; i++)
		if (node_t_check_coaptative_with(nodes[i], bar, nets.count)) field_n[count++] = nodes[i];
	vrtx_t field = {count, field_n};
	return field;
}

vrtx_t get_pow_coapt_field(nets_t nets, int pow){
	int nets_cnt = nets.count;
	data_t data = data_t_construct(nets_cnt,  nets.nets[0].vrtx.count / 10);
	int field_cnt = 0;
	for (int i = 0; i < nets_cnt; i++)
		for (unsigned int j = 0; j < nets.nets[i].vrtx.count; j++){
			node_t* node = nets.nets[i].vrtx.nodes[j];
			if (node_t_get_coaptative_power(node, nets_cnt) == pow){
				add_elem_to_data_ext(data, j, i, nets.nets[i].vrtx.count / 20 + 3);
				field_cnt++;
			}
		}
	node_t** field_n = (node_t**)calloc(field_cnt, sizeof(node_t*));
	int offset = 0;
	for (int i = 0; i < nets_cnt; i++)
		for (int j = 0; j < data.busy_len[i]; j++)
			field_n[offset++] = nets.nets[i].vrtx.nodes[data.data[i][j]];
	vrtx_t field = {(unsigned int)field_cnt, field_n};
	data_t_destruct(data);
	return field;
}

//может не сработать на двух лепестках
node_t* get_next_free_edge_node(node_t** nodes, int cnt, net_t net, int bar){
	node_t* node = nodes[cnt - 1];
	int neigh_cnt = node[0].cnt_springs, i;
	for (i = 0; i < neigh_cnt; i++){
		spring_t* spr = net.springs.springs[node[0].springs_id[i]];
		node_t* neigh = spring_t_get_other_end(spr, node);
		if (!is_free_edge(neigh[0].state) || (nodes[i][0].coaptative > 0 && nodes[i][0].coaptative != bar))
			continue;
		int flag = 0;
		int j;
		for (j = 0; j < cnt; j++)
			if (neigh == nodes[j]) flag = 1;
		if (!flag) return neigh;
	}
	return NULL;
}

vrtx_t net_t_get_free_edge(nets_t nets, int cur, int bar){
	net_t net = nets.nets[cur];
	node_t** nodes = net.vrtx.nodes;
	unsigned int n_cnt = net.vrtx.count, n_edge = 0, i;
	for(i = 0; i < n_cnt; i++)
		if (is_free_edge((*(nodes[i])).state) &&  node_t_check_coaptative_with(nodes[i], bar, nets.count)) n_edge++;
	node_t** free_edge = (node_t**)calloc(n_edge, sizeof(node_t*));
	n_edge = 0;
	for(i = 0; i < n_cnt; i++)
		if (is_free_edge((*(nodes[i])).state) &&  node_t_check_coaptative_with(nodes[i], bar, nets.count))
			free_edge[n_edge++] = nodes[i];

	vrtx_t vrtx = {n_edge, free_edge};
	return vrtx;
}

int node_t_get_coaptative_cnt(net_t net, node_t* node,  int bar, int nets_cnt){
	int neigh_cnt = node[0].cnt_springs, j;
	int coapt_neigh = 0;
	for (j = 0; j < neigh_cnt; j++){
		spring_t* spr = net.springs.springs[node[0].springs_id[j]];
		node_t* neigh = spring_t_get_other_end(spr, node);
		if (node_t_check_coaptative_with(neigh, bar, nets_cnt))
			coapt_neigh++;
	}
	return coapt_neigh;
}

int net_t_set_coapt_state(net_t net, int bar, int nets_cnt){
	int cnt = net.vrtx.count, i;
	int cnt_vrtx = 0;
	for (i = 0; i < cnt; i++){
		node_t* node = net.vrtx.nodes[i];
		node[0].coapt_state = 0;
		if (node_t_check_coaptative_with(node, bar, nets_cnt)){
			node[0].coapt_state = node_t_get_coaptative_cnt(net, node, bar, nets_cnt);
			if (node[0].coapt_state < (int)node[0].cnt_springs) cnt_vrtx++;
		}
	}
	return cnt_vrtx;
}

vrtx_t net_t_get_bottom_contact_bnd(nets_t nets, int cur, int bar){
	net_t net = nets.nets[cur];
	int cnt_vrtx = net_t_set_coapt_state(net, bar, nets.count);
	node_t** bottom_bnd = (node_t**)calloc(cnt_vrtx, sizeof(node_t*));
	int cnt = net.vrtx.count, i, j = 0;
	for (i = 0; i < cnt; i++){
		node_t* node = net.vrtx.nodes[i];
		if (node[0].coapt_state < (int)node[0].cnt_springs && node_t_check_coaptative_with(node, bar, nets.count)) bottom_bnd[j++] = node;
	}

	vrtx_t vrtx = {(unsigned int)j, bottom_bnd};
	return vrtx;
}

double get_distance(node_t* origin, vrtx_t bottom_contact){
	double dist = 1000;
	int id = 0, cnt = bottom_contact.count, i;
	for (i = 0; i < cnt; i++){
		node_t* cur = bottom_contact.nodes[i];
		double cur_dist = SQR_LEN(DIF(cur[0].coord, origin[0].coord));

		if (cur_dist < dist) {
			dist = cur_dist;
			id = i;
		}
	}
	return LEN(DIF(bottom_contact.nodes[id][0].coord, origin[0].coord));
}

int node_in_field(node_t* node, vrtx_t field){
	int j, cnt = field.count;
	for (j = 0; j < cnt; j++)
		if (node == field.nodes[j]) return 1;
	return 0;
}

int elem_in_field(elem_t* elem, vrtx_t field){
	node_t** vrts = elem[0].vrts;
	int i;
	int flag = 0;
	for (i = 0; i < 3; i++)
		flag += node_in_field(vrts[i], field);
	return (flag == 3);
}

point_t get_mid_normal(net_t net, vrtx_t field){
	int cnt = field.count, i, j, flag = 0;
	point_t normal = get_zero_point();
	for (i = 0; i < cnt; i++){
		node_t* node = field.nodes[i];
		int e_cnt = node[0].cnt_elems;
		for (j = 0; j < e_cnt; j++){
			elem_t* elem = net.elems.elems[node[0].elems_id[j]];
			if (elem_in_field(elem, field)){
				normal = point_t_sum(normal, elem[0].or_area);
				flag++;
			}
		}
	}
	if (flag)
		normalize(&normal);
	return normal;
}

int check_bnd(net_t net, vrtx_t field, node_t* node){
	int s_cnt = node[0].cnt_springs, j;
	if (is_bnd(node[0].state)) return 1;
	for (j = 0; j < s_cnt; j++){
		spring_t* spr = net.springs.springs[node[0].springs_id[j]];
		if (!node_in_field(spring_t_get_other_end(spr, node), field))
			return 1;
	}
	return 0;
}

vrtx_t get_field_bnd(net_t net, vrtx_t field){
	int cnt = field.count;
	vrtx_t bnd = vrtx_t_construct(cnt);
	int j = 0;
	for (int i = 0; i < cnt; i++){
		node_t* node = field.nodes[i];
		if (check_bnd(net, field, node)) bnd.nodes[j++] = node;
	}
	bnd.count = j;

	return bnd;
}

double get_directed_dist(node_t* node, vrtx_t bound, point_t direction, point_t plate_norm){
	int cnt = bound.count;
	point_t from = node[0].coord;
	point_t nearest[2] = {};
	double r[2] = {2, 3};

	for (int i = 0; i < cnt; i++){
		node_t* cur = bound.nodes[i];
		if (cur == node) continue;
		point_t or_dir = NORM(DIF(cur[0].coord, from));
		double cur_r = point_t_scal_mul(or_dir, direction);
		if (r[0] > cur_r) {
			r[1] = r[0], nearest[1] = nearest[0];
			r[0] = cur_r, nearest[0] = cur[0].coord;
		}
		else if (r[1] > cur_r)
			r[1] = cur_r, nearest[1] = cur[0].coord;
	}

	point_t line = NORM(DIF(nearest[1], nearest[0]));
	double coef = DOT(DIF(from, nearest[0]), plate_norm) / DOT(line, plate_norm);
	point_t res = nearest[0];
	if (coef >= 1)
		res = nearest[1];
	else if (coef > 0) {
		SCAL_S(coef, &line);
		res = SUM(nearest[0], line);
	}

	double full_dist = LEN(DIF(res, from));
	return full_dist;
}

double net_to_net_get_coapt_directed_depth(nets_t nets, int cur, int bar, point_t direction){
	if (point_t_equal(direction, get_zero_point(), DBL_EPSILON)) return -1;
	vrtx_t field = get_coapt_field(nets, cur, bar);
	point_t normal = get_mid_normal(nets.nets[cur], field);
	point_t plate_norm = NORM(CROSS(direction, normal));
	point_t loc_direction = NORM(DIF(direction, DIR_PROJ(direction, normal)));

	vrtx_t bound = get_field_bnd(nets.nets[cur], field);
	int i, cnt = bound.count;
	double dist = 0;
	for (i = 0; i < cnt; i++){
		node_t* node = bound.nodes[i];
		double cur_dist = get_directed_dist(node, bound, loc_direction, plate_norm);
		if (cur_dist > dist) dist = cur_dist;
	}

	free(field.nodes);
	free(bound.nodes);

	return dist;
}

double net_to_net_get_coapt_unorient_depth(nets_t nets, int cur, int bar){
	vrtx_t bottom_contact = net_t_get_bottom_contact_bnd(nets, cur, bar);
	vrtx_t upper_contact = net_t_get_free_edge(nets, cur, bar);
	double depth = 0;
	int up_cnt = upper_contact.count, i;
	for (i = 0; i < up_cnt; i++){
		node_t* origin = upper_contact.nodes[i];
		double distance = get_distance(origin, bottom_contact);
		if (distance > depth) depth = distance;
	}
	free(bottom_contact.nodes);
	free(upper_contact.nodes);

	return depth;
}

double net_to_net_get_coapt_depth(nets_t nets, int cur, int bar, point_t* direction){
	if (direction == NULL)
		return net_to_net_get_coapt_unorient_depth(nets, cur, bar);

	return net_to_net_get_coapt_directed_depth(nets, cur, bar, *direction);
}

double nets_t_get_coapt_depth(nets_t nets, point_t* direction){
	nets_t_set_coaptative(nets);
	unsigned int nets_cnt = nets.count;
	double depth = 0;
	if (nets_cnt < 2) return depth;
	unsigned int bar = 0, orig = 0;
	for (unsigned int i = 0; i < nets_cnt; i++)
		for (unsigned int j = 0; j < nets_cnt; j++){
			if (j == i) continue;
			double cur_depth = net_to_net_get_coapt_depth(nets, i, j, direction);

			if (depth < cur_depth){
				depth = cur_depth;
				orig = i, bar = j;
			}
			printf("%d %d : cur_depth = %lg\n", i, j, cur_depth);
		}
	double neigh_depth = net_to_net_get_coapt_depth(nets, bar, orig, direction);
	//depth /= nets_cnt * (nets_cnt - 1);
	depth = (depth + neigh_depth) / 2;
	return depth;
}

double nets_t_get_coapt_intersect_depth(nets_t nets, int* pow){
	nets_t_set_coaptative(nets);
	vrtx_t field = get_pow_coapt_field(nets, nets.count - 1);
	printf("Intersect power is %d\n", field.count);
	if (pow) *pow = field.count;
	double sqr_dist = 0;
	for (unsigned int i = 0; i < field.count; i++)
		for (unsigned int j = i + 1; j < field.count; j++){
			double sqr_dist_i_j =point_t_sqr_len(point_t_dif(field.nodes[i][0].coord, field.nodes[j][0].coord));
			if (sqr_dist_i_j > sqr_dist) sqr_dist = sqr_dist_i_j;
		}

	free(field.nodes);

	return sqrt(sqr_dist);
}
//#############################################################################
double net_t_get_coapt_area(net_t net){
	unsigned int cnt_elems = net.elems.count, i;
	double area = 0;
	for (i = 0; i < cnt_elems; i++)
		//area += elem_t_get_coapt_area(net.elems.elems[i]);
		if (elem_t_check_coapt(net.elems.elems[i]))
			area += point_t_length((*(net.elems.elems[i])).or_area);
	return area;
}

double nets_t_get_coapt_area_from_to(nets_t nets, int orig, int bar){
	unsigned int cnt_elems = nets.nets[orig].elems.count;
	double area = 0;
	for (unsigned int i = 0; i < cnt_elems; i++){
		elem_t* elem = nets.nets[orig].elems.elems[i];
		if (elem_t_check_coapt_with(elem, bar, nets.count))
			area += point_t_length(elem[0].or_area);
		}
	return area;
}

double nets_t_get_coapt_area_between(nets_t nets, int net_id1, int net_id2){
	double area1 = nets_t_get_coapt_area_from_to(nets, net_id1, net_id2);
	double area2 = nets_t_get_coapt_area_from_to(nets, net_id2, net_id1);
	double area = (area1 + area2) / 2;
	return area;
}

double net_t_get_full_area(net_t net){
	unsigned int cnt_elems = net.elems.count, i;
	double area = 0;
	for (i = 0; i < cnt_elems; i++)
		area += point_t_length((*(net.elems.elems[i])).or_area);
	return area;
}

