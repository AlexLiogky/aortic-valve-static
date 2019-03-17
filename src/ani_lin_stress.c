#include <math.h>
#include <assert.h>
#include "ani_lin_stress.h"

#define __L_1 106000
#define __L_2 569000
#define __L_3 1200000

const double _L_K[3] = {__L_1, __L_2, __L_3};//{137000, 568000, 968000};
const double _L_X[2] = {0.15, 0.3};
const double _T_K[3] = {__L_1, __L_2, __L_3};//{63000, 570000, 1400000};//{116600, 370000, 600000};
const double _T_X[2] = {0.15, 0.3};

double _f_stiff(double m1, double m2, double sqr_cos_phi){
	return sqrt((m1 * m1 - m2 * m2) * sqr_cos_phi + m2 * m2);
}

double _three_linear_function(const double k[3], const double x[2], double arg){
	double res = 0;
	if (arg < x[0])
		res += k[0] * arg;
	else {
		res += k[0] * x[0];
		if (arg < x[1])
			res += k[1] * (arg - x[0]);
		else {
			res += k[1] * (x[1] - x[0]);
			res += k[2] * (arg - x[1]);
		}
	}
	return res;
}

double _get_longitunial_stiff(double strain){
	int sign = (strain > 0) ? 1 : -1;
	double res = sign * _three_linear_function(_L_K, _L_X, fabs(strain));

	return res;
}

double _get_transverse_stiff(double strain){
	int sign = (strain > 0) ? 1 : -1;
	double res = sign * _three_linear_function(_T_K, _T_X, fabs(strain));

	return res;
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

