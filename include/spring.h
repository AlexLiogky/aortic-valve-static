#ifndef _SPR_
#define _SPR_ 

#include "point.h"
#include "node.h"
#include "stress.h"

typedef struct spring_t{
	unsigned int id;
	double l_0; //initial length
	double l;   //current length
	node_t* ends[2]; //nodes -- ends of spring
	point_t direction;	
	double force;
	stress_t stress_params; //for computation spring constant of stiffness
} spring_t;

typedef struct springs_t{
	unsigned int count;
	spring_t** springs;
}springs_t;

//create spring between two nodes, id - unique identifier of spring
spring_t* spring_t_construct(node_t* node1, node_t* node2, unsigned int id);
void spring_t_destruct(spring_t* spr); //destructor

//create dynamically springs container with size = "count" 
springs_t springs_t_construct(int count);
void springs_t_destruct(springs_t obj); //destructor

//check whether node is one of ends of spring
int spring_t_node_belong(spring_t* spring, node_t* node);

//return not "node" end of spring "spr"
node_t* spring_t_get_other_end(spring_t* spr, node_t* node);

//compute deformation of spring
double spring_t_get_deformation(spring_t* spring);

//update deformational force
void spring_t_update_force(spring_t* spring);

//update direction correspondingly to new coords of nodes
void springs_t_update_direction(spring_t* obj);

//update direction and force for all springs in container
void update_springs(springs_t springs);

//debug print of spring
void spring_t_dump(spring_t* obj);

//return useful size of spring for serialization in bytes
int get_spring_size(spring_t* spr);

//return useful size of spring container including content for serialization in bytes
int get_springs_size(springs_t sprs);

//save seralized spring into buffer 
//ATTENTION: buffer must have enough size
//return current seek in buffer
char* serialize_spring(spring_t* spring, char* buffer);

//save seralized spring container with content into buffer 
//ATTENTION: buffer must have enough size
//return current seek in buffer
char* serialize_springs(springs_t sprs, char* buffer);

//return dynamical deserialized spring
//change offset to next object in buffer
//ATTENTION: deserialization without corresponding nodes container is impossible 
//"vrtx" - nodes container got from buffer 
spring_t* deserialize_spring(char* buffer, int* offset, vrtx_t vrtx);

//return dynamical deserialized node container with content
//change offset to next object in buffer
//ATTENTION: deserialization without corresponding nodes container is impossible 
//"vrtx" - nodes container got from buffer 
springs_t deserialize_springs(char* buffer, int *offset, vrtx_t vrtx);

#endif
