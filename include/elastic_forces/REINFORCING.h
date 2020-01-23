#ifndef _ELASTIC_FORCES_REINFORCING_
#define _ELASTIC_FORCES_REINFORCING_

#include "nets.h"

#ifdef __cplusplus
extern "C" {
#endif

int set_elastic_force_REINFORCING(net_t net, double delta, void* prms);
int precomp_elastic_force_REINFORCING(net_t net, double* prms, void** data);

#ifdef __cplusplus
}
#endif

#endif