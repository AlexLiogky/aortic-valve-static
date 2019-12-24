#include "node.h"
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <stdio.h>

node_t* node_t_construct(double x, double y, double z, double thickness, unsigned int state, unsigned int id){
	node_t* obj = (node_t*) calloc(1, sizeof(node_t));
	point_t current = point_t_get_point(x, y, z);
	point_t_cpy_points(&current, &(*obj).coord);
	point_t_cpy_points(&current, &(*obj).next);
	(*obj).h = thickness;
	(*obj).state = state;
	(*obj).coaptative = -1;
	(*obj).contact_elem_id = NULL;
	(*obj).id = id;
	(*obj).coapt_relation = -1;
	point_t_cpy_points(&current, &(*obj).initial);
	return obj;
}

int is_free(unsigned int state){
	return (state == IN || state == MOB_BND);
}

int is_bnd(unsigned int state){
	return (state == FIX_BND || state == MOB_BND || state == EXTR_BND);
}

int is_fix(unsigned int state){
	return (state == FIX_BND || state == EXTR_BND);
}

int is_free_edge(unsigned int state){
	return (state == MOB_BND || state == EXTR_BND);
}

#define FREE_ARR(X)	if (X != NULL) free(X)

void node_t_destruct(node_t* node){
	FREE_ARR((*node).elems_id);
	FREE_ARR((*node).springs_id);
	FREE_ARR((*node).contact_elem_id);
	free(node);
}
#undef FREE_ARR

vrtx_t vrtx_t_construct(int count){
	vrtx_t vrtx = {};
	vrtx.count = count;
	vrtx.nodes = (node_t**) calloc(count, sizeof(node_t*));
	return vrtx;
}

void vrtx_t_destruct(vrtx_t obj){
	int cnt = obj.count, i;
	for (i = 0; i < cnt; i++) node_t_destruct(obj.nodes[i]);
	free(obj.nodes);
}

int get_node_size(node_t* node){
	int node_size = sizeof(double) + sizeof(point_t) \
		+ sizeof(unsigned int) * (4 + node[0].cnt_elems + node[0].cnt_springs);
	return node_size;
}

int get_vrtx_size(vrtx_t vrtx){
	int res = sizeof(unsigned int);
	for (unsigned int i = 0; i < vrtx.count; i++)
		res += get_node_size(vrtx.nodes[i]);
	return res;
}

#define SAVE_DATA(TYPE, NAME)						\
memcpy(&buffer[off], &(NAME), sizeof(TYPE)); 		\
off += sizeof(TYPE)

//id, coord, h, state, cnt_elems, elems_id, cnt_springs, springs_id
char* serialize_node(node_t* node, char* buffer){
	int off = 0;
	SAVE_DATA(unsigned int, node[0].id);
	SAVE_DATA(point_t, node[0].coord);
	SAVE_DATA(double, node[0].h);
	SAVE_DATA(unsigned int, node[0].state);

	SAVE_DATA(unsigned int, node[0].cnt_elems);
	int cnt = node[0].cnt_elems;
	memcpy(&buffer[off], node[0].elems_id, sizeof(unsigned int) * cnt);
	off += sizeof(unsigned int) * cnt;

	SAVE_DATA(unsigned int, node[0].cnt_springs);
	cnt = node[0].cnt_springs;
	memcpy(&buffer[off], node[0].springs_id, sizeof(unsigned int) * cnt);
	off += sizeof(unsigned int) * cnt;

	return buffer + off;
}

char* serialize_vrtx(vrtx_t vrtx, char* buffer){
	int off = 0;
	SAVE_DATA(unsigned int, vrtx.count);
	int cnt = vrtx.count;
	buffer += off;
	for (int i = 0; i < cnt; i++)
		buffer = serialize_node(vrtx.nodes[i], buffer);
	return buffer;
}

#undef SAVE_DATA

#define DOWNLOAD_DATA(TYPE, NAME)								\
TYPE NAME = *((TYPE*) &buffer[*offset]);						\
*offset += sizeof(TYPE);

node_t* deserialize_node(char* buffer, int* offset){
	DOWNLOAD_DATA(unsigned int, id);
	DOWNLOAD_DATA(point_t, coord);
	DOWNLOAD_DATA(double, h);
	DOWNLOAD_DATA(unsigned int, state);
	DOWNLOAD_DATA(unsigned int, cnt_elems);
	unsigned int* elems_id = (unsigned int*)malloc(cnt_elems * sizeof(unsigned int));
	unsigned int i;
	for (i = 0; i < cnt_elems; i++){
		DOWNLOAD_DATA(unsigned int, data);
		elems_id[i] = data;
	}
	DOWNLOAD_DATA(unsigned int, cnt_springs);
	unsigned int* springs_id = (unsigned int*)malloc(cnt_springs * sizeof(unsigned int));
	for (i = 0; i < cnt_springs; i++){
		DOWNLOAD_DATA(unsigned int, data);
		springs_id[i] = data;
	}
	node_t* node = node_t_construct(coord.coord[0], coord.coord[1], coord.coord[2], h, state, id);
	node[0].cnt_elems = cnt_elems;
	node[0].elems_id = elems_id;
	node[0].cnt_springs = cnt_springs;
	node[0].springs_id =springs_id;

	return node;
}

vrtx_t deserialize_vrtx(char* buffer, int* offset){
	DOWNLOAD_DATA(unsigned int, vrtx_cnt);
	vrtx_t vrtx = vrtx_t_construct(vrtx_cnt);
	for (unsigned int i = 0; i < vrtx_cnt; i++)
		vrtx.nodes[i] = deserialize_node(buffer, offset);

	return vrtx;
}

#undef DOWNLOAD_DATA

int check_uint_in(unsigned int* arr, int count, unsigned int num){
	int i;
	for (i = 0; i < count; i++)
		if (arr[i] == num) return 1;
	return 0;
}

void u_int_arr_dump(unsigned int* arr, unsigned int len){
	printf("arr_len = %u\n", len);
	printf("u_int_arr = %p\n", arr);
	if (arr == NULL) return;
	printf("u_int( ");
	unsigned int del = 20, i, j, ln = len / 20 + 1, cnt = 0;
	for (j = 0; j < ln; j++){
		for (i = 0; i < del; i++)
			if (cnt < len) printf("%u ", arr[cnt++]);
		printf("\n");
	}
	printf(")\n");
}

void node_t_dump(node_t* obj){
	printf("node = %p\n", obj);
	if (obj == NULL) return;
	printf("id: %d\n", (*obj).id);
	printf("coord: ");
	point_t_dump((*obj).coord);
	printf("next: ");
	point_t_dump((*obj).next);
	printf("h = %lg\nstate = %u\n", (*obj).h, (*obj).state);
	printf("elems:{\n");
	u_int_arr_dump((*obj).elems_id, (*obj).cnt_elems);
	printf("}\n");
	printf("springs:{\n");
	u_int_arr_dump((*obj).springs_id, (*obj).cnt_springs);
	printf("}\n");
	return;
}

void vrtx_t_dump(vrtx_t* obj){
	printf("vrtx = %p\n", obj);
	if (obj == NULL) return;
	printf("count = %u\n", (*obj).count);
	printf("nodes:{\n");
	unsigned int i;
	for (i = 0; i < (*obj).count; i++){
		node_t_dump((*obj).nodes[i]);
	}
	printf("}\n");
	return;
}

void node_t_set_shared_elems(node_t* node, unsigned int* elems_id, unsigned int cnt_elems){
	assert(elems_id);
	(*node).elems_id = elems_id;
	(*node).cnt_elems = cnt_elems;
	return;
}

void node_t_set_shared_springs(node_t* node, unsigned int* springs_id, unsigned int cnt_springs){
	assert(springs_id);
	(*node).springs_id = springs_id;
	(*node).cnt_springs = cnt_springs;
	return;
}

void update_nodes(vrtx_t vrtx){
	int cnt_vrt = vrtx.count;
	for (int i = 0; i < cnt_vrt; i++){
		point_t_cpy_points(&(*(vrtx.nodes[i])).next, &(*(vrtx.nodes[i])).coord);
	}
}

