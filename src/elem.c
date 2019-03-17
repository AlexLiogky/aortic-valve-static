#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <assert.h>
#include "elem.h"

elem_t* elem_t_construct(node_t* node1, node_t* node2, node_t* node3, unsigned int id){
	elem_t* obj = (elem_t*) calloc(1, sizeof(elem_t));
	(*obj).coef = 1;
	(*obj).vrts[0] = node1;
	(*obj).vrts[1] = node2;
	(*obj).vrts[2] = node3;
	(*obj).id = id;
	(*obj).neighbours_id = NULL;
	(*obj).cnt_neighbours = 0;
	point_t or_area = point_t_or_area((*(*obj).vrts[0]).coord, (*(*obj).vrts[1]).coord, (*(*obj).vrts[2]).coord);
	(*obj).area_0 = point_t_length(or_area);
	point_t_cpy_points(&or_area, &(*obj).or_area);
	point_t_coef_mul((double)(*obj).coef, &(*obj).or_area);
	return obj;
}

elems_t elems_t_construct(int count){
	elems_t elems = {};
	elems.count = count;
	elems.elems = (elem_t**) calloc(count, sizeof(elem_t*));
	return elems;
}


#define FREE_ARR(X)	if (X != NULL) free(X)

void elem_t_destruct(elem_t* elem){
	FREE_ARR((*elem).neighbours_id);
	free(elem);
}

#undef FREE_ARR

void elems_t_destruct(elems_t obj){
	int cnt = obj.count, i;
	for (i = 0; i < cnt; i++) elem_t_destruct(obj.elems[i]);
	free(obj.elems);
}

void elem_t_dump(elem_t* obj){
	printf("elem = %p\n", obj);
	if (obj == NULL) return;
	printf("coef = %d\n", (*obj).coef);
	printf("or_area:{\n");
	point_t_dump((*obj).or_area);
	printf("}\n");
	printf("vrts = %p\n", (*obj).vrts);
	if ((*obj).vrts == NULL) return;
	for (int i = 0; i < 3; i++){
		printf("vrts[%d] = (", i);
		for (int j = 0; j < 3; j++){
			printf("%lg", (*obj).vrts[i][0].coord.coord[j]);
			if (j < 3 - 1) printf(", ");
		}
		printf(") \n");
	}
	printf("\n");
}

void elem_t_get_springs(elem_t* elem, node_t* springs[3][2]){
	springs[0][0] = (*elem).vrts[0];
	springs[0][1] = (*elem).vrts[1];
	springs[1][0] = (*elem).vrts[0];
	springs[1][1] = (*elem).vrts[2];
	springs[2][0] = (*elem).vrts[1];
	springs[2][1] = (*elem).vrts[2];
	return;
}

void elem_t_update_or_area(elem_t* obj){
	assert(obj);
	assert((*obj).vrts);
	assert((*obj).vrts[0] && (*obj).vrts[1] && (*obj).vrts[2]);
	point_t or_area = point_t_or_area((*(*obj).vrts[0]).coord, (*(*obj).vrts[1]).coord, (*(*obj).vrts[2]).coord);
	point_t_cpy_points(&or_area, &(*obj).or_area);
	point_t_coef_mul((double)(*obj).coef, &(*obj).or_area);
}

void elem_t_update_cntr_mass(elem_t* obj){
	assert(obj);
	assert((*obj).vrts);
	assert((*obj).vrts[0] && (*obj).vrts[1] && (*obj).vrts[2]);
	point_t point1 = (*(*obj).vrts[0]).coord;
	point_t point2 = (*(*obj).vrts[1]).coord;
	point_t point3 = (*(*obj).vrts[2]).coord;
	point_t point = point_t_sum(point1, point2);
	point = point_t_sum(point, point3);
	double r = 1;
	point_t_coef_mul(r / 3, &point);
	point_t_cpy_points(&point, &(*obj).cntr_mass);
}

void update_elems(elems_t elems){
	int cnt_elems = elems.count;
	for (int i = 0; i < cnt_elems; i++){
		elem_t* elem = elems.elems[i];
		elem_t_update_or_area(elem);
		elem_t_update_cntr_mass(elem);
	}
}

point_t point_to_elem_projection(point_t point, elem_t* elem, double* sqr_distance){
	int j;
	point_t triangle[3];
	for (j = 0; j < 3; j++)
		triangle[j] = elem[0].vrts[j][0].coord;
	point_t projection = point_to_triangle_projection(point, triangle, sqr_distance);
	return projection;
}

int get_elem_size(elem_t* elem){
	int elem_size = sizeof(unsigned int) * (5 + elem[0].cnt_neighbours) \
		+ sizeof(point_t) * 2 + sizeof(int) + sizeof(double);
	return elem_size;
}

int get_elems_size(elems_t elems){
	int res = sizeof(unsigned int);
	for (unsigned int i = 0; i < elems.count; i++)
		res += get_elem_size(elems.elems[i]);
	return res;
}

#define SAVE_DATA(TYPE, NAME)						\
memcpy(&buffer[off], &(NAME), sizeof(TYPE)); 		\
off += sizeof(TYPE)

//id, vrts, coef, or_area, cntr_mass, cnt_neighbours, neighbours_id
char* serialize_elem(elem_t* elem, char* buffer){
	int off = 0;
	SAVE_DATA(unsigned int, elem[0].id);
	int i;
	for (i = 0; i < 3; i++){
		SAVE_DATA(unsigned int, elem[0].vrts[i][0].id);
	}
	SAVE_DATA(int, elem[0].coef);
	SAVE_DATA(double, elem[0].area_0);
	SAVE_DATA(point_t, elem[0].or_area);
	SAVE_DATA(point_t, elem[0].cntr_mass);
	SAVE_DATA(unsigned int, elem[0].cnt_neighbours);
	int cnt = elem[0].cnt_neighbours;
	memcpy(&buffer[off], elem[0].neighbours_id, sizeof(unsigned int) * cnt);
	off += sizeof(unsigned int) * cnt;
	
	return buffer + off;
}

char* serialize_elems(elems_t elems, char* buffer){
	int off = 0;
	SAVE_DATA(unsigned int, elems.count);
	int cnt = elems.count;
	buffer += off;
	for (int i = 0; i < cnt; i++)
		buffer = serialize_elem(elems.elems[i], buffer);
	return buffer;
}
#undef SAVE_DATA

//десериализация без сети невозможна
#define DOWNLOAD_DATA(TYPE, NAME)								\
TYPE NAME = *((TYPE*) &buffer[*offset]);						\
*offset += sizeof(TYPE);	

elem_t* deserialize_elem(char* buffer, int* offset, vrtx_t vrtx){
	DOWNLOAD_DATA(unsigned int, id);
	node_t* vrts[3];
	unsigned int i;
	for (i = 0; i < 3; i++){
		DOWNLOAD_DATA(unsigned int, data);
		vrts[i] = vrtx.nodes[data];
	}
	DOWNLOAD_DATA(int, coef);
	DOWNLOAD_DATA(double, area_0);
	DOWNLOAD_DATA(point_t, or_area);
	DOWNLOAD_DATA(point_t, cntr_mass);
	DOWNLOAD_DATA(unsigned int, cnt_neighbours);
	unsigned int* neighbours_id = (unsigned int*)malloc(cnt_neighbours * sizeof(unsigned int));
	for (i = 0; i < cnt_neighbours; i++){
		DOWNLOAD_DATA(unsigned int, data);
		neighbours_id[i] = data;
	}
	elem_t* elem = elem_t_construct(vrts[0], vrts[1], vrts[2], id);
	//elem[0].id = id;
	elem[0].coef = coef;
	elem[0].area_0 = area_0;
	elem[0].or_area = or_area;
	elem[0].cntr_mass = cntr_mass;
	elem[0].cnt_neighbours = cnt_neighbours;
	elem[0].neighbours_id = neighbours_id;
	
	return elem;
}

elems_t deserialize_elems(char* buffer, int *offset, vrtx_t vrtx){
	DOWNLOAD_DATA(unsigned int, elem_cnt);
	elems_t elems = elems_t_construct(elem_cnt);
	for (unsigned int i = 0; i < elem_cnt; i++)
		elems.elems[i] = deserialize_elem(buffer, offset, vrtx);
	
	return elems;
}

int line_to_elem_intersection(point_t line[2], elem_t* elem, point_t* intersect){
	point_t triangle[3];
	for (int i = 0; i < 3; i++)
		triangle[i] = elem[0].vrts[i][0].coord;
	return line_to_triangle_intersection(line, triangle, intersect);
}

int elem_to_triangle(elem_t elem, point_t triangle[3]){
	if (!elem.vrts) return 0;
	for (int i = 0; i < 3; i++)
		if (elem.vrts[i])
			triangle[i] = elem.vrts[i][0].coord;
		else return 0;
	
	return 1;
}

#undef DOWNLOAD_DATA

