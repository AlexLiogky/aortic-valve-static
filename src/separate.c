#include <stdio.h>
#include <math.h>
#include "precomputation.h"
#include "separate.h"

int get_state(int state1, int state2){
	if (state1 == state2) return state1;
	if (is_fix(state1) && is_fix(state2))
		return FIX_BND;
	return IN;
}

unsigned int get_new_vrt_id(node_t* node1, node_t* node2, unsigned int shift){
	unsigned int cnt1 = (*node1).cnt_springs, cnt2 = (*node2).cnt_springs, i = 0, j = 0;
	while (i < cnt1 && j < cnt2 && (*node1).springs_id[i] != (*node2).springs_id[j]){
		unsigned int id1 = (*node1).springs_id[i];
		unsigned int id2 = (*node2).springs_id[j];
		if (id1 < id2) i++;
		else j++;
	}
	return (*node1).springs_id[i] + shift;
}

#define ELEM_CNSTR(X, Y, Z)												\
elems.elems[e_off] = elem_t_construct(X, Y, Z, e_off);				\
if (point_t_scal_mul(or_area, elems.elems[e_off][0].or_area) < 0)	\
	elems.elems[e_off][0].coef = -1;									\
elems.elems[e_off][0].area_0 = area_0;								\
e_off++		

#define E_SPR_CNSTR(X, Y, Z)											\
springs.springs[s_off] = spring_t_construct(X, Y, s_off);			\
springs.springs[s_off][0].stress_params = stress_t_similar_cpy(0.5, net.springs.springs[Z][0].stress_params);\
springs.springs[s_off++][0].l_0 = 0.5 * net.springs.springs[Z][0].l_0

#define S_SPR_CNSTR(X)													\
springs.springs[s_off] = spring_t_construct(X, node, s_off);			\
springs.springs[s_off][0].stress_params = stress_t_similar_cpy(0.5, pastspr[0].stress_params);\
springs.springs[s_off++][0].l_0 = l_0;

net_t get_next_hierarchical_net(net_t net){
	unsigned int e_cnt = net.elems.count;
	unsigned int s_cnt = net.springs.count;
	unsigned int v_cnt = net.vrtx.count;
	unsigned int nv = v_cnt + s_cnt;
	unsigned int ns = 2 * s_cnt + 3 * e_cnt;
	unsigned int ne = 4 * e_cnt, i;
	unsigned int v_off = 0, e_off = 0, s_off = 0;
	vrtx_t vrtx = vrtx_t_construct(nv);
	for (i = 0; i < v_cnt; i++){
		node_t node = net.vrtx.nodes[i][0];
		point_t point = node.coord;
		unsigned int state = node.state;
		double h = node.h;
		vrtx.nodes[v_off] = node_t_construct(point.coord[0], point.coord[1], point.coord[2], h, state, v_off);
		v_off++;
	}
	for (i = v_cnt; i < nv; i++){
		spring_t spr = net.springs.springs[i - v_cnt][0];
		point_t point = point_t_sum(spr.ends[0][0].coord, spr.ends[1][0].coord);
		point_t_coef_mul(0.5, &point);
		unsigned int state = get_state(spr.ends[0][0].state, spr.ends[1][0].state);
		double h = 0.5 * (spr.ends[0][0].h + spr.ends[1][0].h);
		vrtx.nodes[v_off] = node_t_construct(point.coord[0], point.coord[1], point.coord[2], h, state, v_off);
		v_off++;
	}
	elems_t elems = elems_t_construct(ne);
	springs_t springs = springs_t_construct(ns); //нужно дописать этот кусок
	for (i = 0; i < e_cnt; i++){
		double area_0 = net.elems.elems[i][0].area_0 / 4;
		point_t or_area = net.elems.elems[i][0].or_area;
		node_t** pastnodes = net.elems.elems[i][0].vrts;
		node_t* nodes[3];
		nodes[0] = vrtx.nodes[pastnodes[0][0].id];
		nodes[1] = vrtx.nodes[pastnodes[1][0].id];
		nodes[2] = vrtx.nodes[pastnodes[2][0].id];
		unsigned int id1 = get_new_vrt_id(pastnodes[0], pastnodes[1], 0);
		unsigned int id2 = get_new_vrt_id(pastnodes[0], pastnodes[2], 0);
		unsigned int id3 = get_new_vrt_id(pastnodes[1], pastnodes[2], 0);
		node_t* node1 = vrtx.nodes[id1 + v_cnt];
		node_t* node2 = vrtx.nodes[id2 + v_cnt];
		node_t* node3 = vrtx.nodes[id3 + v_cnt];
		ELEM_CNSTR(node1, node2, nodes[0]);
		ELEM_CNSTR(node1, node3, nodes[1]);
		ELEM_CNSTR(node2, node3, nodes[2]);
		ELEM_CNSTR(node1, node2, node3);
		E_SPR_CNSTR(node1, node2, id3);
		E_SPR_CNSTR(node2, node3, id1);
		E_SPR_CNSTR(node1, node3, id2);
	}
	for (i = 0; i < s_cnt; i++){
		node_t** pastends = net.springs.springs[i][0].ends;
		node_t* ends[2];
		ends[0] = vrtx.nodes[pastends[0][0].id];
		ends[1] = vrtx.nodes[pastends[1][0].id];
		unsigned int id = i;
		node_t* node = vrtx.nodes[id + v_cnt];
		spring_t* pastspr = net.springs.springs[id]; 
		double l_0 = 0.5 * pastspr[0].l_0;
		S_SPR_CNSTR(ends[0]);
		S_SPR_CNSTR(ends[1]);
	}
	
	net_t newnet = {vrtx, elems, springs};
	return newnet;	
}
#undef ELEM_CNSTR
#undef E_SPR_CNSTR
#undef S_SPR_CNSTR


nets_t create_next_hierarchical_nets(nets_t nets){
	unsigned int count_nets = nets.count, i;
	nets_t newnets = nets_t_get_net(count_nets);
	for (i = 0; i < count_nets; i++){
		newnets.nets[i] = get_next_hierarchical_net(nets.nets[i]);
	}
	
	init_nets(newnets);
	nets_t_recognize_bnds(newnets);
	nets_t_construct_nodes_contact(newnets);	//order is unimportant
	nets_t_set_elems_neighbours(newnets);		//order is unimportant
	set_relax_consts(newnets, 0.02, 0.05);
	
	return newnets;
}

/*nets_t create_next_hierarchical_dyn_nets(nets_t nets){
	unsigned int count_nets = nets.count, i;
	nets_t newnets = nets_t_get_net(count_nets);
	for (i = 0; i < count_nets; i++){
		if (net_is_static(nets.nets[i])) newnets.nets[i] = cpy_net(nets.nets[i]);
		else newnets.nets[i] = get_next_hierarchical_net(nets.nets[i]);
	}
	
	init_nets(newnets);
	nets_t_recognize_bnds(newnets);
	nets_t_construct_nodes_contact(newnets);	//order is unimportant
	nets_t_set_elems_neighbours(newnets);		//order is unimportant
	set_relax_consts(newnets, 0.02, 0.05);
	
	return newnets;
}*/

//############separation################################################
void set_state(unsigned int* state){
	if (*state == 0) *state = 2;
	if (*state == 2) *state = 3;
	if (*state == 3) *state = 2;
}

unsigned int get_sep_state(spring_t spr){
	unsigned int state1 = spr.ends[0][0].state, state2 = spr.ends[1][0].state;
	set_state(&state1), set_state(&state2);
	if (state1 == state2) return state1 ;
	return 2;
}

void save_elem(FILE* fp, int id1, int id2, int id3, point_t normal){
	fprintf(fp, "%d %d %d %.6f %.6f %.6f\n", id1, id2, id3,\
	 normal.coord[0], normal.coord[1], normal.coord[2]);
}

void sep_net(FILE* fp, net_t net){
	unsigned int e_cnt = net.elems.count;
	unsigned int s_cnt = net.springs.count;
	unsigned int v_cnt = net.vrtx.count;
	unsigned int nv = v_cnt + s_cnt;
	unsigned int ne = 4 * e_cnt, i;
	fprintf(fp, "%u %u\n", nv, ne);
	for (i = 0; i < v_cnt; i++){
		node_t node = net.vrtx.nodes[i][0];
		point_t point = node.coord;
		unsigned int state = node.state;
		set_state(&state);
		fprintf(fp, "%.6f %.6f %.6f %u\n", point.coord[0], point.coord[1], point.coord[2], state);
	}
	for (i = 0; i < s_cnt; i++){
		spring_t spr = net.springs.springs[i][0];
		point_t point = point_t_sum(spr.ends[0][0].coord, spr.ends[1][0].coord);
		point_t_coef_mul(0.5, &point);
		unsigned int state = get_sep_state(spr);
		fprintf(fp, "%.6f %.6f %.6f %u\n", point.coord[0], point.coord[1], point.coord[2], state);
	}
	
	for (i = 0; i < e_cnt; i++){
	point_t or_area = net.elems.elems[i][0].or_area;
	node_t** nodes = net.elems.elems[i][0].vrts; 
	unsigned int id1 = get_new_vrt_id(nodes[0], nodes[1], v_cnt) + 1;
	unsigned int id2 = get_new_vrt_id(nodes[0], nodes[2], v_cnt) + 1;
	unsigned int id3 = get_new_vrt_id(nodes[1], nodes[2], v_cnt) + 1;
	save_elem(fp, id1, id2, nodes[0][0].id + 1, or_area);
	save_elem(fp, id1, id3, nodes[1][0].id + 1, or_area);
	save_elem(fp, id2, id3, nodes[2][0].id + 1, or_area);
	save_elem(fp, id1, id2, id3, or_area);
	}
}

void save_separation(nets_t nets, char* f_name){
	FILE* fp = fopen (f_name, "w");
	
	fprintf(fp, "%u\n", nets.count);
	unsigned int i;
	for (i = 0; i < nets.count; i++)
		sep_net(fp, nets.nets[i]);
	
	fclose(fp);
	return;
}
//############separation################################################

