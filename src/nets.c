#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <stdio.h>
#include <assert.h>
#include "point.h"
#include "node.h"
#include "elem.h"
#include "spring.h"
#include "nets.h"

_net_t_data* _net_t_data_construct(){
    return (_net_t_data*)calloc(1, sizeof(_net_t_data));
}

void _net_t_data_destruct(_net_t_data* nd){
    if (nd && nd->data) free(nd->data);
#ifdef INDIVIDUAL_ELASTIC_INFO
    if (nd && nd->elastic_info) {
        elast_info_t_destruct(nd->elastic_info);
        free(nd->elastic_info);
    }
#endif
    free(nd);
}

net_t net_t_get(vrtx_t vrtx, elems_t elems, springs_t springs){
    net_t net = {vrtx, elems, springs, NETS_DYNAMIC, _net_t_data_construct()};

	return net;
}

int net_is_static(net_t net){
	return (net.state == NETS_STATIC);
}

void net_t_set_state(net_t* net, int state){
	net[0].state = state;
}

nets_t nets_t_get_net(int count){
	nets_t nets;
	nets.count = count;
	nets.nets = (net_t*) calloc(count, sizeof(net_t));
	return nets;
}

#ifdef INDIVIDUAL_ELASTIC_INFO
void elast_info_t_destruct(elast_info_t* data){
    if (data && data->precomp)
        free(data->precomp);
}
#endif

void net_t_destruct(net_t net){
	vrtx_t_destruct(net.vrtx);
	elems_t_destruct(net.elems);
	springs_t_destruct(net.springs);
    _net_t_data_destruct(net.nd);
}

void nets_t_destruct(nets_t nets){
	for (unsigned int i = 0; i < nets.count; i++)
		net_t_destruct(nets.nets[i]);
	free(nets.nets);
}

void nets_t_surfacial_free(nets_t nets){
    free(nets.nets);
}

void nets_t_dump(nets_t nets){
	printf("nets {\n");
	printf("count of nets = %d ", nets.count);
	printf("&nets = %p\n", nets.nets);
	printf("}\n");
}

void net_t_set_node_shared_elems(net_t* net){
	int cnt = (*net).vrtx.count, i;

	int e_cnt = (*net).elems.count, j;
	unsigned int* elem_cnt = (unsigned int*) calloc(cnt, sizeof(int));
	for (i = 0; i < cnt; i++)
		for (j = 0; j < e_cnt; j++){
			node_t** vrts = (*net).elems.elems[j][0].vrts;
			node_t* chk = (*net).vrtx.nodes[i];
			if (vrts[0] == chk || vrts[1] == chk || vrts[2] == chk)
				elem_cnt[i]++;
			}
	for (i = 0; i < cnt; i++){
		unsigned int* elem_id = (unsigned int*) calloc(elem_cnt[i], sizeof(int));
		node_t* chk = (*net).vrtx.nodes[i];
		int k = 0;
		for (j = 0; j < e_cnt; j++){
			node_t** vrts = (*net).elems.elems[j][0].vrts;
			if (vrts[0] == chk || vrts[1] == chk || vrts[2] == chk)
				elem_id[k++] = j;
			}
		node_t_set_shared_elems(chk, elem_id, elem_cnt[i]);
	}
	free(elem_cnt);
}

int spr_t_eq(node_t* ends1[2], node_t* ends2[2]){
	if ((ends1[0] == ends2[0] && ends1[1] == ends2[1]) || (ends1[1] == ends2[0] && ends1[0] == ends2[1]))
		return 1;
	return 0;
}

int check_spring_in_springs(node_t* ends[2], springs_t springs, int count){
	int i;
	for (i = 0; i < count; i++){
		if  (spr_t_eq(ends, (*(springs.springs[i])).ends)) return 1;
	}
	return 0;
}

void net_t_set_springs(net_t* net){		//можно ускорить, применяя data_t для хранения имеющихся spring у заданного node
	int e_cnt = (*net).elems.count, i;
	unsigned int s_cnt = 0, k;
	for(i = 0; i < e_cnt; i++){
		node_t* springs[3][2];
		elem_t_get_springs((*net).elems.elems[i], springs);
		for (k = 0; k < 3; k++){
				if (!check_spring_in_springs(springs[k], (*net).springs, s_cnt)){
					(*net).springs.springs[s_cnt] = spring_t_construct(springs[k][0], springs[k][1], s_cnt);
					s_cnt++;
				}
		}
		if (i % 5000 == 0 && i > 0) printf("Setted springs on %d edges\n", i);
	}
	(*net).springs.count = s_cnt;
}

void net_t_set_node_shared_springs(net_t* net){
	int cnt = (*net).vrtx.count, i, j;
	int s_cnt =(*net).springs.count;
	unsigned int* spr_cnt = (unsigned int*) calloc(cnt, sizeof(int));
	for (i = 0; i < cnt; i++)
		for (j = 0; j < s_cnt; j++){
			node_t** ends = (*((*net).springs.springs[j])).ends;
			node_t* chk = (*net).vrtx.nodes[i];
			if(ends[0] == chk || ends[1] == chk)
				spr_cnt[i]++;
			}

	for (i = 0; i < cnt; i++){
		unsigned int* spring_id = (unsigned int*) calloc(spr_cnt[i], sizeof(int));
		int k = 0;
		node_t* chk = (*net).vrtx.nodes[i];
		for (j = 0; j < s_cnt; j++){
			node_t** ends = (*((*net).springs.springs[j])).ends;
			if(ends[0] == chk || ends[1] == chk)
				spring_id[k++] = j;
			}
		node_t_set_shared_springs(chk, spring_id, spr_cnt[i]);
	}
	free(spr_cnt);
}

void update_net(net_t net){
	update_nodes(net.vrtx);
	update_elems(net.elems);
	update_springs(net.springs);
}

void init_net(net_t* net){
	net_t_set_node_shared_elems(net);
	net_t_set_node_shared_springs(net);
	update_net(*net);
}

void init_nets(nets_t nets){
	unsigned int cnt = nets.count, j;
	for (j = 0; j < cnt; j++)
		init_net(&nets.nets[j]);
}

void update_nets(nets_t nets){
	unsigned int cnt = nets.count, j;
	for (j = 0; j < cnt; j++)
		if (!net_is_static(nets.nets[j])) update_net(nets.nets[j]);
}

point_t node_to_elem_projection(node_t* node, elem_t* elem, double* sqr_distance){
	return point_to_elem_projection(node[0].coord, elem, sqr_distance);
}

point_t point_to_net_projection(point_t point, net_t net, double* sqr_distance){
	int e_cnt = net.elems.count, i;
	double sqr_dist = 1e+20;
	point_t projection;
	for (i = 0; i < e_cnt; i++){
		elem_t* elem = net.elems.elems[i];
		double cur_dist = 2e+20;
		point_t cur_proj = point_to_elem_projection(point, elem, &cur_dist);
		if (sqr_dist > cur_dist) {
			sqr_dist = cur_dist;
			projection = cur_proj;
		}
	}
	if (sqr_distance != NULL) *sqr_distance = sqr_dist;
	return projection;
}

point_t node_to_net_projection(node_t* node, net_t net, double* sqr_distance){
	return point_to_net_projection(node[0].coord, net, sqr_distance);
}

void nets_t_construct_nodes_contact(nets_t nets){
	unsigned int i, cnt = nets.count;
	if (cnt < 2) return;
	for (i = 0; i < cnt; i++){
		if (net_is_static(nets.nets[i])) continue;
		node_t** nodes = nets.nets[i].vrtx.nodes;
		unsigned int n_cnt = nets.nets[i].vrtx.count, j;
		for (j = 0; j < n_cnt; j++){
			node_t* node = nodes[j];
			(*node).contact_elem_id = (int*) malloc(cnt * sizeof(int));
			unsigned int k;
			for (k = 0; k < cnt; k++)
				(*node).contact_elem_id[k] = -1;
		}
	}
}

int get_net_size(net_t net){
	int res = sizeof(int);
	res += get_vrtx_size(net.vrtx);
	res += get_elems_size(net.elems);
	res += get_springs_size(net.springs);
	return res;
}

int get_nets_size(nets_t nets){
	int res = sizeof(unsigned int);
	for (unsigned int i = 0; i < nets.count; i++)
		res += get_net_size(nets.nets[i]);
	return res;
}


char* serialize_net(net_t net, char* buffer){
	buffer = serialize_vrtx(net.vrtx, buffer);
	buffer = serialize_elems(net.elems, buffer);
	buffer = serialize_springs(net.springs, buffer);
	memcpy(buffer, &net.state, sizeof(int));
	buffer += sizeof(int);

	return buffer;
}

char* serialize_nets(nets_t nets, int* fullsize){
	int size = get_nets_size(nets);
	if (fullsize != NULL) *fullsize = size;
	char* buf = (char*)calloc(size, sizeof(char*));
	char* buffer = buf;
	unsigned int cnt = nets.count;
	memcpy(buffer, &cnt, sizeof(unsigned int));
	buffer += sizeof(unsigned int);
	for (unsigned int i = 0; i < cnt; i++)
		buffer = serialize_net(nets.nets[i], buffer);

	return buf;
}

#define DOWNLOAD_DATA(TYPE, NAME)								\
TYPE NAME = *((TYPE*) &buffer[*offset]);						\
*offset += sizeof(TYPE);

net_t deserialize_net(char* buffer, int* offset){
	vrtx_t vrtx = deserialize_vrtx(buffer, offset);
	elems_t elems = deserialize_elems(buffer, offset, vrtx);
	springs_t sprs = deserialize_springs(buffer, offset, vrtx);
	DOWNLOAD_DATA(int, state);
	net_t net = net_t_get(vrtx, elems, sprs);
	net_t_set_state(&net, state);

	return net;
}

nets_t deserialize_nets(char* buffer){
	if (buffer == NULL) return nets_t_get_net(0);
	int offset = 0;
	unsigned int count_nets = *((unsigned int*) &buffer[offset]);
	offset += sizeof(unsigned int);
	nets_t nets = nets_t_get_net(count_nets);
	for (unsigned int i = 0; i < count_nets; i++)
		nets.nets[i] = deserialize_net(buffer, &offset);

	nets_t_construct_nodes_contact(nets);

	return nets;
}

#undef DOWNLOAD_DATA

net_t cpy_net(net_t net){
	int len = get_net_size(net);
	char* buffer = (char*)calloc(len, sizeof(char));
	serialize_net(net, buffer);
	int offset = 0;
	net_t copy = deserialize_net(buffer, &offset);
	free(buffer);
	return copy;
}

nets_t cpy_nets(nets_t nets){
	char* buffer = serialize_nets(nets, NULL);
	nets_t copy = deserialize_nets(buffer);
	free(buffer);
	return copy;
}

spring_t* get_shared_spring(net_t net, node_t* n1, node_t* n2)
{
    node_t* n[2] = {n1, n2};
    int ch = 0;
    if (n1->cnt_springs > n2->cnt_springs) ch = 1;
    for (unsigned int i = 0; i < n[ch]->cnt_springs; ++i)
        if (spring_t_node_belong(net.springs.springs[n[ch]->springs_id[i]], n[(ch + 1) % 2]))
            return net.springs.springs[n[ch]->springs_id[i]];
    return NULL;
}

void net_t_set_thickness(net_t net, double thickness){
    for (int i = 0; i < net.vrtx.count; ++i)
        net.vrtx.nodes[i]->h = thickness;
}

