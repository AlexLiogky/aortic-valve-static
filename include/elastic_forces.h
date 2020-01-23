#ifndef _ELASTIC_FORCES_
#define _ELASTIC_FORCES_

#include <stdlib.h>
#include "elastic_forces/MSM.h"
#include "elastic_forces/TRQS.h"
#include "elastic_forces/REINFORCING.h"
#include "elastic_forces/NEOGOOK.h"

#ifdef __cplusplus
extern "C" {
#endif

double* resize(void** data, long nmembs);
double* get_params(void* data);
int dummy_saving_prms(void** data, double* prms, int cnt);
int set_elastic_force_dummy(net_t net, double delta, void* prms, point_t (*f)(net_t net, node_t* node, double* prms));

#ifdef __cplusplus
}
#endif


#endif