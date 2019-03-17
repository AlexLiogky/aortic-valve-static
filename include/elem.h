#ifndef _ELEM_
#define _ELEM_ 

#include "point.h"
#include "node.h"

typedef struct elem_t{
	unsigned int id;
	node_t* vrts[3]; //nodes -- vertices
	int coef; //set to -1 if compute internal normal
	double area_0;
	point_t or_area;
	point_t cntr_mass;
	unsigned int cnt_neighbours;
	unsigned int* neighbours_id;
} elem_t;

typedef struct elems_t{
	unsigned int count;
	elem_t** elems;
}elems_t;

//create dynamicaly triangle element on 3 nodes, id - unique identifier of elem
elem_t* elem_t_construct(node_t* node1, node_t* node2, node_t* node3, unsigned int id);
void elem_t_destruct(elem_t* elem); //destructor

//create dynamicaly elems container with size = "count" 
elems_t elems_t_construct(int count);
void elems_t_destruct(elems_t obj); //destructor

//return boundary edges in "springs[3][2]"
void elem_t_get_springs(elem_t* elem, node_t* springs[3][2]);

//debug print of elem
void elem_t_dump(elem_t* obj);

//update or_area correspondingly to nodes coords
void elem_t_update_or_area(elem_t* obj);

//update cntr_mass correspondingly to nodes coords
void elem_t_update_cntr_mass(elem_t* obj);

//update or_area and cntr_mass for all elems in container
void update_elems(elems_t elems);

//return projection of point to a triangle surface of element
//if (sqr_dist != NULL) save a square of distance from point to projection here
point_t point_to_elem_projection(point_t point, elem_t* elem, double* sqr_distance);

//return useful size of elem for serialization in bytes
int get_elem_size(elem_t* elem);

//return useful size of elem container including content for serialization in bytes
int get_elems_size(elems_t elems);

//save seralized elem into buffer 
//ATTENTION: buffer must have enough size
//return current seek in buffer_id
char* serialize_elem(elem_t* elem, char* buffer);

//save seralized elem container with content into buffer 
//ATTENTION: buffer must have enough size
//return current seek in buffer
char* serialize_elems(elems_t elems, char* buffer);

//return dynamical deserialized elem
//change offset to next object in buffer
//ATTENTION: deserialization without corresponding nodes container is impossible 
//"vrtx" - nodes container got from buffe
elem_t* deserialize_elem(char* buffer, int* offset, vrtx_t vrtx);

//return dynamical deserialized elem container with content
//change offset to next object in buffer
//ATTENTION: deserialization without corresponding nodes container is impossible 
//"vrtx" - nodes container got from buffer 
elems_t deserialize_elems(char* buffer, int *offset, vrtx_t vrtx);

//save intersection between line and triangle surface of element into "intersect"
//return count of intersections (2 means infinite count of intersections)
int line_to_elem_intersection(point_t line[2], elem_t* elem, point_t* intersect);

//convert elem to triangle and save result into "triangle[3]"
//return 1 on success
int elem_to_triangle(elem_t elem, point_t triangle[3]);

#endif
