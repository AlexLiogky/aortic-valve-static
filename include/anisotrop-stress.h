#ifndef ANISOTROP_STRESS_H_INCLUDED
#define ANISOTROP_STRESS_H_INCLUDED

#include "stress-interface.h"

typedef struct stress_t{
	double sigma1;
	double sigma2;
	double sigma3;
} stress_t;

#endif // ANISOTROP-STRESS_H_INCLUDED
