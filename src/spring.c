#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "spring.h"

spring_t* spring_t_construct(node_t* node1, node_t* node2, unsigned int id){
	spring_t* obj = (spring_t*) calloc(1, sizeof(spring_t));
	(*obj).ends[0] = node1;
	(*obj).ends[1] = node2;
	point_t direction = point_t_dif((*node2).coord, (*node1).coord);
	point_t_cpy_points(&direction, &(*obj).direction);
	double len = point_t_length(direction);
	(*obj).l_0 = len;
	(*obj).l = len;
	point_t_coef_mul(1 / len, &(*obj).direction);
	(*obj).force = 0;
	(*obj).id = id;
	return obj;
}

springs_t springs_t_construct(int count){
	springs_t springs = {};
	springs.count = count;
	springs.springs = (spring_t**) calloc(count, sizeof(spring_t*));
	return springs;
}

void spring_t_destruct(spring_t* spr){
	free(spr);
}

void springs_t_destruct(springs_t obj){
	int cnt = obj.count, i;
	for (i = 0; i < cnt; i++) spring_t_destruct(obj.springs[i]);
	free(obj.springs);
}

void spring_t_dump(spring_t* obj){
	printf("spring = %p\n", obj);
	if (obj == NULL) return;
	printf("direction: ");
	point_t_dump((*obj).direction);
	printf("l_0 = %lg\n", (*obj).l_0);
	printf("l = %lg\n", (*obj).l);
	printf("ends = %p", (*obj).ends);
	if ((*obj).ends == NULL) {
		printf("\n");
		return;
	}
	printf("{\n");
	int i;
	for (i = 0; i < 2; i++)
		node_t_dump((*obj).ends[i]);
	printf("}\n");
}

int spring_t_node_belong(spring_t* spring, node_t* node){
	 return ((*spring).ends[0] == node || (*spring).ends[1] == node);
}

node_t* spring_t_get_other_end(spring_t* spr, node_t* node){
	return (node != (*spr).ends[0]) ? (*spr).ends[0] : (*spr).ends[1];
}

double spring_t_get_deformation(spring_t* spring){
	double l = (*spring).l, l_0 = (*spring).l_0;
	double eps = (l - l_0) / l_0;
	return eps;
}

void spring_t_update_force(spring_t* spring){
	double eps = spring_t_get_deformation(spring);
	double abs_f = stress_t_get_force(eps, (*spring).stress_params);
	(*spring).force = abs_f;
}

void springs_t_update_direction(spring_t* obj){
	node_t* node1 = (*obj).ends[0];
	node_t* node2 = (*obj).ends[1];
	point_t direction = point_t_dif((*node2).coord, (*node1).coord);
	double len = point_t_length(direction);
	(*obj).l = len;
	point_t_cpy_points(&direction, &(*obj).direction);
	point_t_coef_mul(1 / len, &(*obj).direction);
}

void update_springs(springs_t springs){
	int cnt_sprs = springs.count;
	for (int i = 0; i < cnt_sprs; i++){
		spring_t* spr = springs.springs[i];
		springs_t_update_direction(spr);
		spring_t_update_force(spr);
	}
}

int get_spring_size(spring_t* spr){
	int spr_size = sizeof(unsigned int) * 3 + 3 * sizeof(double) + sizeof(point_t) + sizeof(stress_t);
	return spr_size;
}

int get_springs_size(springs_t sprs){
	int res = sizeof(unsigned int);
	for (unsigned int i = 0; i < sprs.count; i++)
		res += get_spring_size(sprs.springs[i]);
	return res;
}

#define SAVE_DATA(TYPE, NAME)						\
memcpy(&buffer[off], &(NAME), sizeof(TYPE)); 		\
off += sizeof(TYPE)

//id, l_0, l, ends, direction, force, stress_params
char* serialize_spring(spring_t* spring, char* buffer){
	int off = 0;
	SAVE_DATA(unsigned int, spring[0].id);
	SAVE_DATA(double, spring[0].l_0);
	SAVE_DATA(double, spring[0].l);
	for (int i = 0; i < 2; i++){
		SAVE_DATA(unsigned int, spring[0].ends[i][0].id);
	}
	SAVE_DATA(point_t, spring[0].direction);	
	SAVE_DATA(double, spring[0].force);
	SAVE_DATA(stress_t, spring[0].stress_params);
	
	return buffer + off;
}

char* serialize_springs(springs_t sprs, char* buffer){
	int off = 0;
	SAVE_DATA(unsigned int, sprs.count);
	int cnt = sprs.count;
	buffer += off;
	for (int i = 0; i < cnt; i++)
		buffer = serialize_spring(sprs.springs[i], buffer);
	return buffer;
}
#undef SAVE_DATA

//десериализация без сети невозможна
#define DOWNLOAD_DATA(TYPE, NAME)								\
TYPE NAME = *((TYPE*) &buffer[*offset]);						\
*offset += sizeof(TYPE);

spring_t* deserialize_spring(char* buffer, int* offset, vrtx_t vrtx){
	DOWNLOAD_DATA(unsigned int, id);
	DOWNLOAD_DATA(double, l_0);
	DOWNLOAD_DATA(double, l);
	node_t* ends[2];
	int i;
	for (i = 0; i < 2; i++){
		DOWNLOAD_DATA(unsigned int, data);
		ends[i] = vrtx.nodes[data];
	}
	DOWNLOAD_DATA(point_t, direction);
	DOWNLOAD_DATA(double, force);
	DOWNLOAD_DATA(stress_t, stress_params);
	//stress_t stress_params = {0};
	//*offset += 2 * sizeof(double);
	spring_t* spring = spring_t_construct(ends[0], ends[1], id);
	spring[0].l_0 = l_0;
	spring[0].l = l;
	spring[0].direction = direction;
	spring[0].force = force;
	spring[0].stress_params = stress_params;
	
	return spring;
}

springs_t deserialize_springs(char* buffer, int *offset, vrtx_t vrtx){
	DOWNLOAD_DATA(unsigned int, sprs_cnt);
	springs_t sprs = springs_t_construct(sprs_cnt);
	for (unsigned int i = 0; i < sprs_cnt; i++)
		sprs.springs[i] = deserialize_spring(buffer, offset, vrtx);
		
	return sprs;
}

#undef DOWNLOAD_DATA

