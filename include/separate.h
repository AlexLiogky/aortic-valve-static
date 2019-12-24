#ifndef _SEPARATOR_H
#define _SEPARATOR_H

#include "nets.h"

#ifdef __cplusplus
extern "C" {
#endif

//return a new net with half the grid spacing,
//while updating the spring stiffness parameters
net_t get_next_hierarchical_net(net_t net);

//return a new net with half the grid spacing,
//while updating the spring stiffness parameters
nets_t create_next_hierarchical_nets(nets_t nets);


//save nets to file in special format
void save_separation(nets_t nets, char* f_name);
//############separation################################################

#ifdef __cplusplus
}
#endif


#endif
