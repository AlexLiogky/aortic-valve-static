#include <math.h>
#include "anisotrop-stress.h"

const double _E1_0 = 137000;	//Pa
const double _E1_1 = 568000;	//Pa
const double _Eps_lim1 = 0.175;
const double _E2_0 = 63000;	//Pa
const double _E2_1 = 570000;	//Pa
const double _Eps_lim2 = 0.175;

double _f_stiff(double m1, double m2, double sqr_cos_phi){
	return sqrt((m1 * m1 - m2 * m2) * sqr_cos_phi + m2 * m2);
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
	double sigma1 = _f_stiff(_E1_0, _E2_0, cos_phi);
	double sigma2 = (_Eps_lim1 > _Eps_lim2) ? _f_stiff(_E1_0, _E2_1, cos_phi) : _f_stiff(_E1_1, _E2_0, cos_phi);
	double sigma3 = _f_stiff(_E1_1, _E2_1, cos_phi);
	sigma1 *= coef;
	sigma2 *= coef;
	sigma3 *= coef;
	stress_t stress = {sigma1, sigma2, sigma3};
	return stress;
}

stress_t stress_t_similar_cpy(double coef, stress_t stress){
	double sigma1 = stress.sigma1 * coef;
	double sigma2 = stress.sigma2 * coef;
	double sigma3 = stress.sigma3 * coef;
	stress_t newstress = {sigma1, sigma2, sigma3};
	return newstress;
}

double stress_t_get_force(double deformation, stress_t stress){
	double force = 0;
	double eps_lim1 = (_Eps_lim1 < _Eps_lim2) ? _Eps_lim1 : _Eps_lim2;
	double eps_lim2 = (_Eps_lim1 < _Eps_lim2) ? _Eps_lim2 : _Eps_lim1;
	double sigma1 = stress.sigma1;
	double sigma2 = stress.sigma2;
	double sigma3 = stress.sigma3;
	int sign = (deformation > 0) ? 1 : -1;
	force += (fabs(deformation) < eps_lim1) ? sigma1 * deformation : sigma1 * eps_lim1 * sign;
	if (fabs(deformation) > eps_lim1){
		double deform = deformation - sign * eps_lim1;
		force += (fabs(deformation) < eps_lim2) ? sigma2 * deform  : sigma2 * eps_lim2 * sign;
	}
	force += (fabs(deformation) > eps_lim2) ? sigma3 * (deformation - sign * eps_lim2) : 0;
	return force;
}

