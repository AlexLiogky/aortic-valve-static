#include <math.h>
#include <assert.h>
#include "isotrop-stress.h"

///жёсткости увеличены в 100 раз

const double _E_0 = 1000000;//10000000;		//Pa в 100
const double _E_1 = 1000000;//70600000;	//Pa в 10
const double _Eps_lim = 0.18;

stress_t stress_t_construct(point_t spr_direction, point_t* fiber_direction, double coef){
	double sigma1 = _E_0;
	double sigma2 = _E_1;
	sigma1 *= coef;
	sigma2 *= coef;
	stress_t stress = {sigma1, sigma2};
	return stress;
}

stress_t stress_t_similar_cpy(double coef, stress_t stress){
	double sigma1 = stress.sigma1 * coef;
	double sigma2 = stress.sigma2 * coef;
	stress_t newstress = {sigma1, sigma2};
	return newstress;
}

double stress_t_get_force(double deformation, stress_t stress){
	double force = 0;
	double eps_lim = _Eps_lim;
	double sigma1 = stress.sigma1;
	double sigma2 = stress.sigma2;
	int sign = (deformation > 0) ? 1 : -1;
	force += (fabs(deformation) < eps_lim) ? sigma1 * deformation : sigma1 * eps_lim * sign;
	force += (fabs(deformation) > eps_lim) ? sigma2 * (deformation - sign * eps_lim) : 0;
	return force;
}

