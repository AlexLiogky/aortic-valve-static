#ifndef _COMPUTE_VALVE_H
#define _COMPUTE_VALVE_H

#include "nets.h"
#include "bound-box.h"
#include <sys/time.h>

#ifdef __cplusplus
extern "C" {
#endif

#define FLAG_REL 1
extern double Contact_Force_Const1;
extern double Contact_Force_Const2;
extern double gt_elastic;

enum ElasticModels{
    EMOD_MSM,
    EMOD_TRQS,
    EMOD_NEOGOOK,
    EMOD_REINFORCING,
    EMOD_COUNT
};

typedef struct def_elast_info_t{
    int model;
    double params[MAX_COUNT_ELASTIC_PARAMS];
    int version;
} def_elast_info_t;

//############computation###############################################
void set_default_elastic_model(int type);
void set_default_elastic_params(double* prms, int count);
#ifdef INDIVIDUAL_ELASTIC_INFO
void set_elastic_info(net_t* net, int type, double* prms, int count);
#endif

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


typedef struct collision_data_t{
    box_t box;
    int check_freq;
    int counter;
} collision_t;

typedef struct solver_data_t{
    double delta;           //coefficient converting force into shift
    double eps;             //stop condition
    int ElasticModelType;   //number of used elastic model
} solver_t;

typedef struct statistical_data_t{
    point_t* mid_shift;
    unsigned int cnt_nodes;
    int it_cnt;
    double init_div;
} statistic_t;

typedef struct initial_world_condition_t{
    double P;
} wrld_cnd_t;

typedef struct world_t{
    nets_t dynamic_nets;
    nets_t static_nets;
    nets_t union_nets;
    collision_t collision;
    wrld_cnd_t conditions;
    solver_t solver_data;
    statistic_t statistic;
} world_t;

world_t* world_t_construct(nets_t dynamic_nets, nets_t static_nets, solver_t solver_data, collision_t collision_data, wrld_cnd_t conditions);

statistic_t statistical_data_construct(nets_t nets);
void statistical_data_destruct(statistic_t *st);

nets_t create_union_net(nets_t nets1, nets_t nets2);

collision_t collision_data_t_construct(nets_t dynamic_nets, nets_t static_nets, int check_freq);
void collision_data_t_destruct(collision_t cl);

void collision_data_t_update(nets_t nets, collision_t* coll_data);

int collision_t_get_state(collision_t coll_data);

void update_statistic_t(nets_t nets, statistic_t *st);

double statistic_t_get_max_mid_diviation(statistic_t st);

double statistic_t_get_full_mid_diviation(statistic_t st);

void statistic_t_reset(statistic_t *st);

void set_initial_solving_params(world_t* world);
long double compute_dif_ms(struct timeval end, struct timeval start);

void compute_nets_time(long double compute_time, world_t* world, int max_its);

void relaxation(nets_t nets, double coef);
//############computation###############################################

#ifdef __cplusplus
}
#endif


#endif
