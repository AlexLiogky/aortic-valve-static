#ifndef ISOTROP_STRESS_H
#define ISOTROP_STRESS_H

#include "stress-interface.h"

#ifdef __cplusplus
extern "C" {
#endif

typedef struct stress_t{
	double sigma1;
	double sigma2;
	//double sigma3;
} stress_t;

#ifdef __cplusplus
}
#endif


#endif
