#include <math.h>
#include "ani_exp_stress.h"

#define __L_A 6000
#define __L_B 11

const double _L_A = __L_A;//9617.7;//1.932;
const double _L_B = __L_B;//10;//51.6;
const double _T_A = __L_A;//1017.3;//9007.2;//0.078;
const double _T_B = __L_B;//16.97;//9.822;//13.2;

double _f_stiff(double m1, double m2, double sqr_cos_phi){
	return sqrt((m1 * m1 - m2 * m2) * sqr_cos_phi + m2 * m2);
}

double _get_longitunial_stiff(double deformation){
	const double a = _L_A, b = _L_B;
	return a * sinh(b * deformation);
}

double _get_transverse_stiff(double deformation){
	const double a = _T_A, b = _T_B;
	return a * sinh(b * deformation);
}

stress_t stress_t_construct(point_t spr_direction, point_t* fiber_direction, double coef){
	assert(fiber_direction);
	point_t fib_dir;
	point_t_cpy_points(fiber_direction, &fib_dir);
	//normalize(&fib_dir);
	//normalize(&spr_direction);
	double cos_phi = DOT(fib_dir, spr_direction);
	double sqr_cos_phi = cos_phi * cos_phi / \
							(SQR_LEN(fib_dir) * SQR_LEN(spr_direction));
	stress_t stress = {sqr_cos_phi, coef};
	return stress;
}

stress_t stress_t_similar_cpy(double coef, stress_t stress){
	double sqr_cos_phi = stress.sqr_cos_phi;
	double newcoef = stress.coef * coef;
	stress_t newstress = {sqr_cos_phi, newcoef};
	return newstress;
}

double stress_t_get_force(double deformation, stress_t stress){
	double m1 = _get_longitunial_stiff(deformation);
	double m2 = _get_transverse_stiff(deformation);
	int sign = (deformation > 0) ? 1 : -1;
	double force = sign * stress.coef * _f_stiff(m1, m2, stress.sqr_cos_phi);

	return force;
}

