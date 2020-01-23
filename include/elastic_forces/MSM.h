#ifndef _ELASTIC_FORCES_MSM_
#define _ELASTIC_FORCES_MSM_

#include "nets.h"

#ifdef __cplusplus
extern "C" {
#endif
int set_elastic_force_MSM(net_t net, double delta, void* prms);
int precomp_elastic_force_MSM(net_t net, double* prms, void** data);

#ifdef __cplusplus
}
#endif

#endif