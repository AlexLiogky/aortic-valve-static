#ifndef _NETS_
#define _NETS_

#include "point.h"
#include "node.h"
#include "elem.h"
#include "spring.h"

#ifdef __cplusplus
extern "C" {
#endif

typedef struct net_t{
	vrtx_t vrtx;
	elems_t elems;
	springs_t springs;
	int state;
}net_t;

typedef struct nets_t{
	net_t* nets;
	unsigned int count;
}nets_t;

enum NetsMoveState{
	NETS_DYNAMIC,
	NETS_STATIC
};

//######################################################################
//######################################################################
//create net from prepared parts
//ATTENTION: only create union object, doesn't set required links between parts
net_t net_t_get(vrtx_t vrtx, elems_t elems, springs_t springs);

//check move state, return true if static
int net_is_static(net_t net);

void net_t_set_state(net_t* net, int state);

//create dynamicaly net container
nets_t nets_t_get_net(int count);

void net_t_destruct(net_t net); //destructor

//container destructor, delete content too
void nets_t_destruct(nets_t nets);

//container destructor, doesn't change content
void nets_t_surfacial_free(nets_t nets);

//debug print of nets container
void nets_t_dump(nets_t nets);

//set shared elems to all nodes in net
void net_t_set_node_shared_elems(net_t* net);

//checks whether the ends of the two springs are equal
int spr_t_eq(node_t* ends1[2], node_t* ends2[2]);

//check whether "ends" belong to any spring in springs container
int check_spring_in_springs(node_t* ends[2], springs_t springs, int count);

//compute all springs produced by elements and set it into "net"
void net_t_set_springs(net_t* net);

//set shared springs to all nodes in net
void net_t_set_node_shared_springs(net_t* net);

//update net correspondingly to new coords of nodes
void update_net(net_t net);

//set neccassary links between parts of net
void init_net(net_t* net);

//init all nets in net container
void init_nets(nets_t nets);

//update all nets in net container correspondingly to new coords of nodes
void update_nets(nets_t nets);

//return projection of node coord to elem surface
//if (sqr_distance != NULL) set square of distance between "node" and "elem" here
point_t node_to_elem_projection(node_t* node, elem_t* elem, double* sqr_distance);

//return the nearest point of net surface to "point"
//if (sqr_distance != NULL) set square of distance between "point" and "net" here
point_t point_to_net_projection(point_t point, net_t net, double* sqr_distance);

//return the nearest point of net surface to "node.coord"
//if (sqr_distance != NULL) set square of distance between "node" and "net" here
point_t node_to_net_projection(node_t* node, net_t net, double* sqr_distance);

//initialize some optimizational collision parameters
void nets_t_construct_nodes_contact(nets_t nets);

//return useful size of net for serialization in bytes
//including serialization size of content
int get_net_size(net_t net);

//return useful size of net container including content for serialization in bytes
int get_nets_size(nets_t nets);

//save seralized net into buffer
//ATTENTION: buffer must have enough size
//return current seek in buffer_id
char* serialize_net(net_t net, char* buffer);

//save seralized net container with content into buffer
//ATTENTION: buffer must have enough size
//return current seek in buffer
char* serialize_nets(nets_t nets, int* fullsize);

//return dynamical deserialized net from buffer
//change "offset" to next object in buffer
net_t deserialize_net(char* buffer, int* offset);

//return dynamical deserialized net container with content from buffer
nets_t deserialize_nets(char* buffer);

//create dynamical deep copy of "net"
net_t cpy_net(net_t net);

//create dynamical deep copy of net container
nets_t cpy_nets(nets_t nets);

//return spring between n1 and n2
//if there is no such spring then return NULL
spring_t* get_shared_spring(net_t net, node_t* n1, node_t* n2);

#ifdef __cplusplus
}
#endif


#endif
