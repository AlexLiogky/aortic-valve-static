#ifndef STRESS_INTERFACE_H_INCLUDED
#define STRESS_INTERFACE_H_INCLUDED

#include "point.h"

typedef struct stress_t stress_tt;

//construct a spring stress params
//ATTENTION: current state is supposed undeformed
stress_tt stress_t_construct(point_t spr_direction, point_t* fiber_direction, double coef);

//copy scaled stress params
//it's useful for creating hierachical net
//coef - scale coefficient, relation length of new spring to the old
stress_tt stress_t_similar_cpy(double coef, stress_tt stress);

//retrun force caused by deformation
double stress_t_get_force(double deformation, stress_tt stress);

#endif // STRESS_INTERFACE_H_INCLUDED
