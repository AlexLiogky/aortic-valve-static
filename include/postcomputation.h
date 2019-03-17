#ifndef _POSTCOMPUTE_VALVE_H
#define _POSTCOMPUTE_VALVE_H

#include "nets.h"

//return middle area of "node"
double node_t_area(net_t net, node_t* node);

//print some statistical information about nets
void print_statistic(nets_t nets);

//return depth of coaptation between "cur" and "bar"
//suppose that depth of coaptation is maximal distance between 
//points of bottom and the least deviating point of upper bounds from "direction" of coaptation field 
double net_to_net_get_coapt_directed_depth(nets_t nets, int cur, int bar, point_t direction);

//return depth of coaptation between "cur" and "bar"
//suppose that depth of coaptation is maximal distance between 
//points of bottom and upper bounds of coaptation field
double net_to_net_get_coapt_unorient_depth(nets_t nets, int cur, int bar);

//return depth of coaptation between "cur" and "bar"
//shape function choosing "directed_depth" or "unorient_depth" 
//depending on "direction" =? NULL 
double net_to_net_get_coapt_depth(nets_t nets, int cur, int bar, point_t* direction);

//return middle depth of coaptation
double nets_t_get_coapt_depth(nets_t nets, point_t* direction);

//return height of central hole in aortic valve
double nets_t_get_coapt_intersect_depth(nets_t nets);

//return total area of coaptation on "net"
double net_t_get_coapt_area(net_t net);

//return area of coaptation from "orig" to "bar"
double nets_t_get_coapt_area_from_to(nets_t nets, int orig, int bar);

//return middle area of coaptation between "orig" to "bar"
double nets_t_get_coapt_area_between(nets_t nets, int net_id1, int net_id2);

//return full area of "net"
double net_t_get_full_area(net_t net);

#endif
