#include <stdlib.h>
#include <math.h>
#include <stdio.h>
#include <string.h>
#include <float.h>
#include "point.h"

point_t* point_t_construct(double x, double y, double z){
	point_t* point = (point_t*) calloc(1, sizeof(point_t));
	(*point).coord[0] = x;
	(*point).coord[1] = y;
	(*point).coord[2] = z;
	return point;
}

point_t point_t_get_point(double x, double y, double z){
	point_t point;
	point.coord[0] = x;
	point.coord[1] = y;
	point.coord[2] = z;
	
	return point;
}

int point_t_equal(point_t point1, point_t point2, double err){
	int i, res = 1;
	for (i = 0; i < DIM; i++){
		res *= (fabs(point1.coord[i] - point2.coord[i]) < err);
	}
	
	return res;
}

void point_t_coef_mul(double coef, point_t* point){
	int i;
	for (i = 0; i < DIM; i++)
		(*point).coord[i] *= coef;
}

point_t point_t_coef_mul_new(double coef, point_t point){
	point_t res;
	for (int i = 0; i < DIM; i++)
		res.coord[i] = point.coord[i] * coef;
	
	return res;
}

point_t point_t_sum(point_t point1, point_t point2){
	point_t sum;
	int i;
	for (i = 0; i < DIM; i++){
		sum.coord[i] = (point1).coord[i] + (point2).coord[i];
	}

	return sum;
}

point_t point_t_add_scale(point_t origin, double coef, point_t addition){
	point_t res;
	for (int i = 0; i < DIM; i++)
		res.coord[i] = origin.coord[i] + coef * addition.coord[i];
	return res;
}

void point_t_add_scale_self(point_t* origin, double coef, point_t addition){
	for (int i = 0; i < DIM; i++)
		origin[0].coord[i] = origin[0].coord[i] + coef * addition.coord[i];
}

point_t point_t_dif(point_t point1, point_t point2){
	point_t dif;
	int i;
	for (i = 0; i < DIM; i++)
		dif.coord[i] = point1.coord[i] - point2.coord[i];
	return dif;
}

double point_t_scal_mul(point_t point1, point_t point2){
	double res = 0;
	int i;
	for (i = 0; i < DIM; i++)
		res += point1.coord[i] * point2.coord[i];
	return res;
}

double point_t_sqr_len(point_t point){
	return point_t_scal_mul(point, point);
}

double point_t_length(point_t point){
	return sqrt(point_t_scal_mul(point, point));
}

point_t point_t_vec_mul(point_t point1, point_t point2){
	double	x1 = point1.coord[0], x2 = point2.coord[0],\
			y1 = point1.coord[1], y2 = point2.coord[1],\
			z1 = point1.coord[2], z2 = point2.coord[2];
	double 	x = y1*z2 - y2*z1,\
			y = x2*z1 - x1*z2,\
			z = x1*y2 - x2*y1;
	return point_t_get_point(x, y, z);
}

void normalize(point_t* point){
	point_t_coef_mul(1.0 / point_t_length(*point), point);
}

point_t normalize_new(point_t point){
	return point_t_coef_mul_new(1.0 / point_t_length(point), point);
}

point_t get_zero_point(){
	return point_t_get_point(0, 0, 0);
}

int point_t_less(point_t pnt1, point_t pnt2){
	for (int i = 0; i < DIM; i++)
		if (pnt1.coord[i] >= pnt2.coord[i]) return 0;
	return 1;
}

int point_t_less_eq(point_t pnt1, point_t pnt2){
	for (int i = 0; i < DIM; i++)
		if (pnt1.coord[i] > pnt2.coord[i]) return 0;
	return 1;
}

int point_t_greater(point_t pnt1, point_t pnt2){
	for (int i = 0; i < DIM; i++)
		if (pnt1.coord[i] <= pnt2.coord[i]) return 0;
	return 1;
}

int point_t_greater_eq(point_t pnt1, point_t pnt2){
	for (int i = 0; i < DIM; i++)
		if (pnt1.coord[i] < pnt2.coord[i]) return 0;
	return 1;
}

point_t point_t_or_area(point_t point1, point_t point2, point_t point3){
	point_t triag_area = SCAL(0.5, CROSS(DIF(point2, point1), DIF(point3, point1)));
	return triag_area;
}

point_t point_t_cpy(point_t point0){
	point_t point;
	int i;
	for (i = 0; i < DIM; i++)
		point.coord[i] = point0.coord[i];
	return point;
}

void point_t_cpy_points(point_t* src, point_t* tgt){
	memcpy((*tgt).coord, (*src).coord, DIM * sizeof(double));
}

double get_area(point_t point1, point_t point2, point_t point3){
	return LEN(OR_AREA(point1, point2, point3));
}

void point_t_dump(point_t point){ 
	double x = point.coord[0], y = point.coord[1], z = point.coord[2];
	printf("point = (%lg, %lg, %lg)\n", x, y, z);
}

point_t signed_norm_inf(point_t point, int* sign, int* n_max_crd){ //(x, 1, z)
	double max_comp = fabs(point.coord[0]);
	int max_j = 0;
	for (int i = 1; i < DIM; i++){
		double abs_comp = fabs(point.coord[i]);
		if (max_comp < abs_comp) {
			max_comp = abs_comp;
			max_j = i;
		}
	}
	if (sign == NULL || *sign != 0) max_comp = point.coord[max_j];
	if (n_max_crd != NULL) *n_max_crd = max_j;
	return SCAL(1 / max_comp, point);
}

point_t point_to_direction_projection(point_t vec, point_t dir){
	return SCAL(DOT(dir, vec) / SQR_LEN(dir), dir);
}

point_t point_to_line_projection(point_t point, point_t line[2]){
	point_t dir = DIF(line[1], line[0]);
	point_t rel = DIF(point, line[0]);
	point_t projection = SUM(line[0], DIR_PROJ(rel, dir));
	return projection;
}

point_t point_to_line_segment_projection(point_t point, point_t line_segment[2]){
	point_t direction = DIF(line_segment[1], line_segment[0]);
	point_t deflection = DIF(point, line_segment[0]);
	double coef = DOT(direction, deflection);
	if (coef < 0) return line_segment[0];
	coef /= SQR_LEN(direction);
	if (coef > 1) return line_segment[1];
	point_t projection = ADD(line_segment[0], coef, direction);
	return projection;
}

//если хранить нормаль у элементов, то следующую функцию можно 
//оптимизировать для вычисления элементов
point_t point_to_plate_projection(point_t point, point_t plate[3]){
	point_t normal = OR_AREA(plate[0], plate[1], plate[2]);
	point_t proj = SUM(point, DIR_PROJ(DIF(plate[0], point), normal));
	return proj;
}

#define _S(X, Y, Z) a.coord[X]*b.coord[Y]*c.coord[Z]
double point_t_or_volume(point_t point1, point_t point2, point_t point3, point_t point4){
	point_t a = point_t_dif(point2, point1);
	point_t b = point_t_dif(point3, point2);
	point_t c = point_t_dif(point4, point3);
	double volume = _S(0, 1, 2) + _S(1, 2, 0) + _S(2, 0, 1) - _S(2, 1, 0) - _S(1, 0, 2) - _S(0, 2, 1);
	return volume;
}
#undef _S

int _sign(double num){
	return (num > 0) ? 1 : -1;
}

int is_parallel(point_t point1, point_t point2){
	point_t point3 = point_t_vec_mul(point1, point2);
	return EQ_ERR(point3, ZERO(), 5 * DBL_EPSILON);
}

int is_ortogonal(point_t point1, point_t point2){
	return (fabs(point_t_scal_mul(point1, point2)) < 10 * DBL_EPSILON);
}

int line_to_normal_plate_intersection(point_t dir, point_t line, point_t normal, point_t plate, point_t* intersect){
	if (!is_ortogonal(normal, dir)){
		if (intersect != NULL)
			*intersect = ADD(line, DOT(DIF(plate, line), normal) / DOT(dir, normal), dir);
		return 1;
	}else{
		if (fabs(DOT(normal, DIF(plate, line))) > 10 * DBL_EPSILON)
			return 0;
		if (intersect != NULL) *intersect = line;
		return 2;	
	}	
}

int line_to_plate_intersection(point_t line[2], point_t plate[3], point_t* intersect){
	point_t normal = OR_AREA(plate[0], plate[1], plate[2]);
	point_t dir = DIF(line[1], line[0]);
	return line_to_normal_plate_intersection(dir, line[0], normal, plate[0], intersect);
}

int is_point_belongs_triangle(point_t point, point_t triangle[3]){//point must belong to triangle plate
	point_t sides[3][2];
	sides[0][0] = triangle[0], sides[0][1] = triangle[1];
	sides[1][0] = triangle[1], sides[1][1] = triangle[2];
	sides[2][0] = triangle[2], sides[2][1] = triangle[0];
	point_t normals[3];
	for (int i = 0; i < 3; i++)
		normals[i] = CROSS(DIF(sides[i][1], sides[i][0]), DIF(point, sides[i][1]));
	return (DOT(normals[0], normals[1]) > 0 && DOT(normals[1], normals[2]) > 0);
}

int line_intersect_line(point_t l1[2], point_t l2[2], point_t* intersect){ //lines must belong to one plate
	point_t d2 = DIF(l2[1], l2[0]), d1 = DIF(l1[1], l1[0]);
	point_t normal = DIF(d2, DIR_PROJ(d2, d1));
	
	if (!EQ_ERR(normal, ZERO(), 10 * DBL_EPSILON)){
		return line_to_normal_plate_intersection(d2, l2[0], normal, l1[0], intersect);
	}else{
		if (is_parallel(DIF(l2[1], l1[0]), d1)){
			if (intersect != NULL) *intersect = l1[0];
			return 2;
		}
		return 0;
	}
}

int line_intersect_line_segment(point_t l[2], point_t l_seg[2], point_t* intersect){ //lines must belong to one plate
	point_t cur_inter;
	int state = line_intersect_line(l, l_seg, &cur_inter);
	if (state == 0) return 0;
	if (state == 1) {
		point_t shift = DIF(cur_inter, l_seg[0]);
		point_t d = DIF(l_seg[1], l_seg[0]);
		double res = DOT(shift, d) / SQR_LEN(d);
		if (res > 0 && res < 1) {
			if (intersect != NULL) *intersect = cur_inter;
			return 1;
		}
		return 0;
	}
	if (intersect != NULL) *intersect = l_seg[0];
	return 2;
}

int line_to_triangle_intersection(point_t line[2], point_t triangle[3], point_t* intersect){
	point_t cur_inter;
	int state = line_to_plate_intersection(line, triangle, &cur_inter);
	if (state == 0) return 0;
	if (state == 1) {
		if (is_point_belongs_triangle(cur_inter, triangle)){
			*intersect = cur_inter;
			return 1;
		}
		return 0;
	}
	point_t sides[3][2];
	sides[0][0] = triangle[0], sides[0][1] = triangle[1];
	sides[1][0] = triangle[1], sides[1][1] = triangle[2];
	sides[2][0] = triangle[2], sides[2][1] = triangle[0];
	for (int i = 0; i < 3; i++)
		if (line_intersect_line_segment(line, sides[i], intersect) > 0) 
			return 2;
	return 0;
}

//эта функция рассматривает пересечение с внутренностью треугольника
int line_fragment_to_trangle_intersection(point_t start, point_t shift, point_t triangle[4], point_t* intersect){ //triangle[3] = normal
	point_t normal = triangle[3];
	double param1 = DOT(DIF(triangle[0], start), normal);
	double param2 = DOT(shift, normal);
	if (fabs(param2) < DBL_EPSILON) {
		if (intersect && fabs(param1) < DBL_EPSILON) *intersect = start;
		return (fabs(param1) < DBL_EPSILON);
	}
	double coef = param1 / param2;
	point_t cur_inter;
	if (coef <= 1 && coef >= 0){
		cur_inter = ADD(start, coef, shift);
	}
	else return 0;
	if (is_point_belongs_triangle(cur_inter, triangle)){
		if (intersect) *intersect = cur_inter;
		return 1;
	}
	
	return 0;
}

point_t point_to_triangle_projection(point_t point, point_t triangle[3], double* sqr_dist){
	point_t plate_projection = point_to_plate_projection(point, triangle);
	point_t sides[3][2];
	sides[0][0] = triangle[0], sides[0][1] = triangle[1];
	sides[1][0] = triangle[1], sides[1][1] = triangle[2];
	sides[2][0] = triangle[2], sides[2][1] = triangle[0];
	point_t normals[3];
	for (int i = 0; i < 3; i++)
		normals[i] = CROSS(DIF(sides[i][1], sides[i][0]), DIF(plate_projection, sides[i][1]));
	if (DOT(normals[0], normals[1]) > 0 && DOT(normals[1], normals[2]) > 0){
		if (sqr_dist != NULL)
			*sqr_dist = SQR_LEN(DIF(plate_projection, point));
		return plate_projection;
	}
	
	point_t projs[3];
	double min_proj = 1e+20;
	int ip;
	for (int i = 0; i < 3; i++){
		projs[i] = point_to_line_segment_projection(point, sides[i]);
		point_t dif = DIF(point, projs[i]);
		double dist = DOT(dif, dif);
		if (dist < min_proj){
			min_proj = dist;
			ip = i;
		}
	}
	if (sqr_dist != NULL) {
			*sqr_dist = min_proj;
		}
	
	return projs[ip];
	
}



matrix3D_t matrix3D_t_get_from_points(point_t points[DIM]){
	matrix3D_t matr;
	for (int i = 0; i < DIM; i++)
		for (int j = 0; j < DIM; j++)
			matr.matrix[i][j] = points[j].coord[i];
	return matr;
}

matrix3D_t matrix3D_t_mul(matrix3D_t m1, matrix3D_t m2){
	matrix3D_t m3 = {0};
	for (int i = 0; i < DIM; i++)
		for (int j = 0; j < DIM; j++)
			for (int k = 0; k < DIM; k++)
				m3.matrix[i][j] += m1.matrix[i][k] * m2.matrix[k][j];
	return m3;
}

double matrix3D_t_det(matrix3D_t m){
	double det = 0;
	for (int i = 0; i < 3; i++)
		det +=  m.matrix[0][i] * (m.matrix[1][(i + 1) % 3] * m.matrix[2][(i + 2) % 3] -\
				m.matrix[1][(i + 2) % 3] * m.matrix[2][(i + 1) % 3]);
	return det;
}

matrix3D_t matrix3D_t_inverse(matrix3D_t m){
	matrix3D_t res;
	double det = matrix3D_t_det(m);
	for(int i = 0; i < 3; i++)
      for(int j = 0; j < 3; j++)
           res.matrix[j][i] = (m.matrix[(i + 1) % 3][(j + 1) % 3] * m.matrix[(i + 2) % 3][(j + 2) % 3] -\
						m.matrix[(i + 1) % 3][(j + 2) % 3] * m.matrix[(i + 2) % 3][(j + 1) % 3]) / det;
	
	return res;
}

void matrix3D_t_dump(matrix3D_t m){
	for (int i = 0; i < DIM; i++){
		for (int j = 0; j < DIM; j++)
			printf("%.3f ", m.matrix[i][j]);
		printf("\n");
	}
}

point_t matrix3D_t_vec_mul(matrix3D_t mat, point_t vec){
	point_t res = {0};
	for (int i = 0; i < DIM; i++)
		for(int j = 0; j < DIM; j++)
			res.coord[i] += mat.matrix[i][j] * vec.coord[j];
	return res;
}


f_point_t f_point_t_convert(point_t point){
	f_point_t mypoint;
	for (int i = 0; i < DIM; i++)
		mypoint.coord[i] = (float) point.coord[i];
	
	return mypoint;
}

point_t f_point_t_to_point(f_point_t p){
	point_t res;
	for (int i = 0; i < DIM; ++i)
		res.coord[i] = p.coord[i];
	return res;
}

int f_point_t_eq(f_point_t p1, f_point_t p2){
	int flag = 0;
	for (int i = 0; i < 3; i++)
		flag += (fabsf(p1.coord[i] - p2.coord[i]) <= FLT_EPSILON);
	return (flag == 3);
} 

void line_t_dump(line_t line){
	printf("line_t:: { len = %d\n", line.pnt_cnt);
	for (int i = 0; i < line.pnt_cnt; i++)
		{printf("%d: ", i); point_t_dump(line.line[i]);}
	printf("}\n");
}

