#ifndef _PRECOMPUTE_VALVE_
#define _PRECOMPUTE_VALVE_

#include "nets.h"
#include <math.h>

#ifdef __cplusplus
extern "C" {
#endif

extern double Contact_Resolution;
extern double Contact_Const;
extern double Min_Contact_Const;

extern double Allow_shift;
extern double Max_shift;

extern point_t Fiber_Direction;
extern double** __Dihedral_Angles;
//############precomputation############################################
//return max value
double max_3d(double x1, double x2, double x3);

//set some collision constants based on nets geometry
void auto_set_contact_recognition_consts(nets_t nets, double* min_resolution);

//set bounds state to bounded nodes of net
void net_t_recognize_bnds(net_t net);

//set bounds state to bounded nodes of nets
void nets_t_recognize_bnds(nets_t nets);

//set elastic parameters to springs as for intial relaxed state
void net_t_set_relax_state(net_t net);

//set elastic parameters to springs as for intial relaxed state
//to all net in nets, direction of anisotropy "fiber_dir" suppose common
void nets_t_set_relax_state(nets_t nets, point_t fiber_dir);

//set shared elems for every elem in net
void net_t_set_elems_neighbours(net_t net);
//set shared elems only for elem in net, advanced low-level function, using doesn't advised
void net_t_elem_t_set_elems_neighbours(net_t net, elem_t* elem);

//set shared elem for every elem in net container
void nets_t_set_elems_neighbours(nets_t nets);

//set some relaxational consts bounded max shift during computations
void set_relax_consts(nets_t nets, double rel_allow_shift, double rel_max_shift);

void set_contact_recognition_resolution(double contact_resolution);
void set_contact_to_plate_recognition_const(double contact_const);
void set_min_contact_recognition_const(double min_contact_const);
void set_contact_recognition_consts(double contact_resolution, double contact_const, double min_contact_const);

double get_Contact_Resolution();
double get_Contact_Const();
double get_Min_Contact_Const();

//do all nessaccary precomputations to prepare nets for computation
//don't set springs parameters
void precomputation(nets_t nets);
//############precomputation############################################

int net_t_count_shared_elems(node_t* node1, node_t* node2, int* elems_id);

#ifdef __cplusplus
}
#endif


#endif
