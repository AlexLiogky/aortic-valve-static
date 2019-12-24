#ifndef _NODE_
#define _NODE_

#include "point.h"

#ifdef __cplusplus
extern "C" {
#endif

typedef struct node_t{
	unsigned int id;
	point_t coord; //current coordinates
	point_t next;
	double h; //thickness
	unsigned int state; //mobility, boundary?
	unsigned int cnt_elems;
	unsigned int* elems_id; //in ascending order!!!
	unsigned int cnt_springs;
	unsigned int* springs_id; //in ascending order!!!
	int coaptative;
	int coapt_state;
	int* contact_elem_id;
	double coapt_relation;
	point_t initial;
}  node_t;

typedef struct vrtx_t{
	unsigned int count;
	node_t** nodes;
}vrtx_t;

enum State
{
    IN = 0, 		//internal, mobile
    FIX_BND = 1,	//boundary, unmobile
    MOB_BND = 2,	//boundary, mobile
    EXTR_BND = 3	//boundary, unmobile, belongs to free edge
};

//dinamicaly construct the node
// x,y,z - current coords
// thickness - local thickness in area of node
// state - node value of State
// id - unique identificator of node
node_t* node_t_construct(double x, double y, double z, double thickness, unsigned int state, unsigned int id);
void node_t_destruct(node_t* node); // destructor

//check state functions
int is_free(unsigned int state);	//is node free?
int is_bnd(unsigned int state);		//is node boundary?
int is_fix(unsigned int state);		//is node fixed?
int is_free_edge(unsigned int state);//is node belongs to free edge?

//create dynamically nodes container with size = "count"
vrtx_t vrtx_t_construct(int count);
void vrtx_t_destruct(vrtx_t obj); // destructor

//return useful size of node for serialization in bytes
int get_node_size(node_t* node);

//return size of node container including content for serialization in bytes
int get_vrtx_size(vrtx_t vrtx);

//save seralized node into buffer
//ATTENTION: buffer must have enough size
//return current seek in buffer
char* serialize_node(node_t* node, char* buffer);

//save seralized node container with content into buffer
//ATTENTION: buffer must have enough size
//return current seek in buffer
char* serialize_vrtx(vrtx_t vrtx, char* buffer);

//return dynamical deserialized node
//change offset to next object in buffer
node_t* deserialize_node(char* buffer, int* offset);

//return dynamical deserialized node container with content
//change offset to next object in buffer
vrtx_t deserialize_vrtx(char* buffer, int* offset);

//checks whether a "num" belongs to an array "arr" of len = "count"
int check_uint_in(unsigned int* arr, int count, unsigned int num);

//debug print of unsigned int array
void u_int_arr_dump(unsigned int* arr, unsigned int len);

//debug print of node
void node_t_dump(node_t* obj);

//debug print of node container
void vrtx_t_dump(vrtx_t* obj);

//set elems_id as shared elements of node
//ATTENTION: set pointer, there is no copy
void node_t_set_shared_elems(node_t* node, unsigned int* elems_id, unsigned int cnt_elems);

//set springs_id as shared springs of node
//ATTENTION: set pointer, there is no copy
void node_t_set_shared_springs(node_t* node, unsigned int* springs_id, unsigned int cnt_springs);

//place next -> coord
void update_nodes(vrtx_t vrtx);

#ifdef __cplusplus
}
#endif


#endif
