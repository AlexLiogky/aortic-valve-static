#ifndef _COMPUTE_VALVE_H
#define _COMPUTE_VALVE_H

#include "nets.h"
#include "bound-box.h"

#define FLAG_REL 1
extern double Contact_Force_Const1;
extern double Contact_Force_Const2;

//############computation###############################################
//set specific collision force
void set_contact_force_consts(double k1, double k2);

//return computed pressure force acted to one node
point_t pressure_force(net_t net, node_t* node, double P);

//experimental function that fit optimal delta for solver
double recommended_delta(nets_t nets, double P);

//return strain force acted to "node" from "spring"
point_t f_spring(net_t net, node_t* node, spring_t* spring);

//return total strain force acted to "node" from all springs
point_t elastic_force(net_t net, node_t* node);

//return contact force acted to node 
//force - total force acted on node without contact
//normal - external normal of contacting body
//dist - approximated distance between node and contacting body
point_t get_contact_force(double dist, point_t normal, double force, node_t* node);

//check whether "node" and "elem" contact now
//f - total force acted on node without contact
//if (r) set distance between node and element center mass 
int check_contact(elem_t* elem, node_t* node, point_t* f, double* r);

//return contact force acted to "node" belonging to "orig" net in "nets"
//from "bar" net in "nets"
point_t contact_force(nets_t nets, int orig, int bar, node_t* node);

//return contact force acted to "node" belonging to "orig" net in "nets"
//from "bar" net in "nets"
//optimazed version of function "contact_force", use bounded box	
point_t box_contact_force(nets_t nets, int orig, int bar, node_t* node, box_t box);

//compute forces without regard bodies collision
//return max ||x_new - x_old||_{infty}, maximal shift over node
double compute_free_nexts(net_t net, double P, double delta);

double get_diviation(nets_t nets); //return |x_new - x_old|

//compute forces with regard bodies collision
//if (flag != 0) set contact force to zero
double compute_constraint_nexts(nets_t nets, int flag, box_t box);

//compute forces acted to every nodes of all nets with regard bodies collision
//if (constraint != 0) set contact force to zero
//make some action preventing unstable behaviour of method
int compute_nexts(nets_t nets, double P, double delta, int constraint, box_t box);	

//compute n iterations of system pseudoevolution
//P - external pressure, suppose P > 0
//delta - method coefficient, that manage convertation of forces into shifts, shift = delta * force
//eps - the computational accuracy at which calculations are stopped, if |current force| < eps * |initial force|
//freq - after how many iterations to output debugging information 
void compute_nets(nets_t nets, double P, double delta, unsigned int n, double eps, int freq);

//make a linear converting to come from "init_points" to "final_points"
//scale and rotate net
void deform_net(net_t net, point_t init_points[3], point_t final_points[3]);

//return first from "line[0]" to "line[1]" intersection between line and net
//if (id) save here id of intersected element
// max_sqr_cntr_mass_dist - geometrical parameter of net, it's possible get it with "get_max_sqr_cntr_mass_dist"
point_t find_line_to_net_intersection(net_t net, point_t line[2], double max_sqr_cntr_mass_dist, int* id);

//return square of max distance between vertex of element and element's center mass
double get_max_sqr_cntr_mass_dist(net_t net);

//move fixed points of "flat" across normal vector to the intersection with "net"
//some algorithm of sewing the flat object to a net 
void flat_obj_projection_to_net(net_t flat, net_t net);

//return len of free edge of net
double net_t_get_len_free_edge(net_t net);

//return len of fix edge of net
double net_t_get_len_fix_edge(net_t net);

//sew net "leaflet" to line "bnd"
//start_node_id - id of first node which belongs to stitched boundary of "leaflet"
void sew_leaflet_to_bnd(net_t leaflet, unsigned int start_node_id, line_t bnd);

//sew leaflet to aorta
//init_attachment - points of commisure on leaflet
//init direction - direction from center of "init_attachment" to the bottom point of leaflet
//attachment - points of commisure on aorta
//blood_direction - approximate direction of central symmetry axis in aorta
//start_node_id - id of first node which belongs to stitched boundary of "leaflet"
//bnd - line on aorta for sewing
void sew_leaflet_to_aorta(	net_t leaflet, point_t init_attachment[2], point_t init_direction,\
							net_t aorta,   point_t attachment[2], 	   point_t blood_direction,\
							unsigned int start_node_id, line_t bnd);
//############computation###############################################

#endif
