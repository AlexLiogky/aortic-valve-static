#ifndef _POINT_H
#define _POINT_H

#ifdef __cplusplus
extern "C" {
#endif

#define DIM 3

typedef struct point_t{
	double coord[DIM];
}point_t;

point_t* point_t_construct(double x, double y, double z); 	//alloc memory and create a point
point_t point_t_get_point(double x, double y, double z); 	//return local variable
point_t get_zero_point(); 									//return (0,0,0)
point_t point_t_cpy(point_t point0); 						//return copy of the point
void point_t_cpy_points(point_t* src, point_t* tgt); 		//return dynamic allocated copy of the point
void point_t_coef_mul(double coef, point_t* point); 		//x <- a*x
point_t point_t_coef_mul_new(double coef, point_t point); 	//return y = a*x
point_t point_t_sum(point_t point1, point_t point2);		//return z = y + x
point_t point_t_add_scale(point_t origin, double coef, point_t addition); //return z = x + a * y
void point_t_add_scale_self(point_t* origin, double coef, point_t addition); //x <- x + a * y
point_t point_t_dif(point_t point1, point_t point2); 		//return z = x - y
double point_t_scal_mul(point_t point1, point_t point2); 	//return (x, y)
double point_t_sqr_len(point_t point); 						//return |x|^2
double point_t_length(point_t point); 						//return |x|
point_t point_t_vec_mul(point_t point1, point_t point2); 	//return [x, y]
void normalize(point_t* point); 							//x <- x / |x|
point_t normalize_new(point_t point); 						//return y = x / |x|
point_t point_t_scale_add_scale(double coef1, point_t p1, double coef2, point_t p2); //return z = c1 * y + c2 * x

//check the equility of points with accuracy "err"
int point_t_equal(point_t point1, point_t point2, double err);

//return 1 if every coord of x is less than correspondig coord of y else 0
int point_t_less(point_t pnt1, point_t pnt2);

//return 1 if every coord of x is less of equal than correspondig coord of y else 0
int point_t_less_eq(point_t pnt1, point_t pnt2);

//return 1 if every coord of x is greater than correspondig coord of y else 0
int point_t_greater(point_t pnt1, point_t pnt2);

//return 1 if every coord of x is greater or equal than correspondig coord of y else 0
int point_t_greater_eq(point_t pnt1, point_t pnt2);

//return oriented area of triangle built on the points
point_t point_t_or_area(point_t point1, point_t point2, point_t point3);

//return area of triangle built on the points
double get_area(point_t point1, point_t point2, point_t point3);

//return (p2 - p1, p3 - p2, p4 - p3)
double point_t_or_volume(point_t point1, point_t point2, point_t point3, point_t point4);

//debug print to standart output content of the point
void point_t_dump(point_t point);

//retrun x / |x[i_inf]| * sig(x[i_inf]),
//where i_inf = arg max_i |x[i]| and
//		sig(x[i_inf]) = 1,				if !(*sign)
//		sig(x[i_inf]) = sign(x[i_inf]), else
point_t signed_norm_inf(point_t point, int* sign, int* n_max_crd);

//return projection of vector vec to axis dir
point_t point_to_direction_projection(point_t vec, point_t dir);

//return projection of point to the line defined by two points
//the checking of nonequlity of points doesn't produce
point_t point_to_line_projection(point_t point, point_t line[2]);

//return projection of point to line sigment defined by their ends
//the checking of nonequlity of points doesn't produce
point_t point_to_line_segment_projection(point_t point, point_t line_segment[2]);

//return projection of point to plate surface defined by three points
//the non-degeneracy of the triangle is not checked
point_t point_to_plate_projection(point_t point, point_t plate[3]);

//return 1 if vectors is parallel
int is_parallel(point_t point1, point_t point2);

//return 1 if vectors is orthogonal
int is_ortogonal(point_t point1, point_t point2);

//save intersection's point between line and plate in intersect
//line defined by own point "line" and directional vector "dir"
//plate defined by own point "plate" and normal vector "normal"
//return count of intersections (2 means line belongs plate)
int line_to_normal_plate_intersection(point_t dir, point_t line, point_t normal, point_t plate, point_t* intersect);

//save intersection's point between line and plate in intersect
//line defined by own three points
//plate defined by own point "plate" and normal vector "normal"
//return count of intersections (2 means line belongs plate)
int line_to_plate_intersection(point_t line[2], point_t plate[3], point_t* intersect);

//checks if a point belongs to a triangle in the affine space of a triangle
//ATTENTION: point must belong to triangle plate
int is_point_belongs_triangle(point_t point, point_t triangle[3]);

//save intersection between lines lying in the same plane into "intersect"
//lines defined by own two points
//return count of intersections (2 means lines coincide)
//ATTENTION: lines must belong to one plate
int line_intersect_line(point_t l1[2], point_t l2[2], point_t* intersect);

//save intersection between line and line segment lying in the same plane into "intersect"
//line defined by own two points
//line segment defined by their ends
//return count of intersections (2 means line segment belongs to line)
//ATTENTION: line and line segment must belong to one plate
int line_intersect_line_segment(point_t l[2], point_t l_seg[2], point_t* intersect);

//save intersection between line and triangle into "intersect"
//return count of intersections (2 means infinite count of intersections)
int line_to_triangle_intersection(point_t line[2], point_t triangle[3], point_t* intersect);

//save intersection between line segment and interior of triangle into "intersect"
//line defined by own point "start" and directional vector "shift"
//triangle defined by own 3 points and for optimization
//"triangle[3]" must contain normal vector
//return count of intersections (2 means infinite count of intersections)
int line_fragment_to_trangle_intersection(point_t start, point_t shift, point_t triangle[4], point_t* intersect);

//return projection of point to a triangle
//if (sqr_dist != NULL) save a square of distance from point to projection here
point_t point_to_triangle_projection(point_t point, point_t triangle[3], double* sqr_dist);

typedef struct matrix3D_t{
	double matrix[DIM][DIM];
}matrix3D_t;

//construct matrix from vectors as columns
matrix3D_t matrix3D_t_get_from_points(point_t points[DIM]); //return |p1 p2 p3|
matrix3D_t matrix3D_t_mul(matrix3D_t m1, matrix3D_t m2); //return A * B
double matrix3D_t_det(matrix3D_t m); //return det(A)
matrix3D_t matrix3D_t_inverse(matrix3D_t m); //return A^-1
point_t matrix3D_t_vec_mul(matrix3D_t mat, point_t vec); //return A * x

//debug print to standart output content of the matrix
void matrix3D_t_dump(matrix3D_t m);

#define DOT(vec1, vec2) point_t_scal_mul(vec1, vec2)		//	(a, b)
#define CROSS(vec1, vec2) point_t_vec_mul(vec1, vec2)		//	[a x b]
#define SCAL(alpha, vec) point_t_coef_mul_new(alpha, vec)	//	k*a
#define SCAL_S(alpha, vec_adr) point_t_coef_mul(alpha, vec_adr)//a <- k*a
#define NORM_S(vec_adr) normalize(vec_adr)					//	a/|a|
#define NORM(vec) normalize_new(vec)						//	a/|a|
#define SQR_LEN(vec) point_t_sqr_len(vec)					//  |a|^2
#define LEN(vec) point_t_length(vec)						//  |a|
#define SUM(vec1, vec2) point_t_sum(vec1, vec2)				//  a + b
#define DIF(vec1, vec2) point_t_dif(vec1, vec2)				//  a - b
#define EQ(vec1, vec2) point_t_equal(vec1, vec2, DBL_EPSILON)// a == b
#define EQ_ERR(vec1, vec2, err) point_t_equal(vec1, vec2, err)//a == b
#define ADD(vec1, alpha, vec2) point_t_add_scale(vec1, alpha, vec2)// a + k*b
#define ADD_S(vec1_adr, alpha, vec2) point_t_add_scale_self(vec1_adr, alpha, vec2)
#define SCAL_SUM(coef1, vec1, coef2, vec2) point_t_scale_add_scale(coef1, vec1, coef2, vec2)
#define VEC(X, Y, Z) point_t_get_point(X, Y, Z)				//(x, y, z)
#define ZERO() get_zero_point()								//(0, 0, 0)
#define LNE(vec1, vec2) point_t_less(vec1, vec2)			//  a < b
#define LEQ(vec1, vec2) point_t_less_eq(vec1, vec2)			// a <= b
#define GNE(vec1, vec2) point_t_greater(vec1, vec2)			// 	a > b
#define GEQ(vec1, vec2) point_t_greater_eq(vec1, vec2)		// a >= b
#define OR_AREA(v1, v2, v3) point_t_or_area(v1, v2, v3)		//[p2 - p1, p3 - p1] / 2
#define DIR_PROJ(vec, dir) point_to_direction_projection(vec, dir)

typedef struct f_point_t{
	float coord[DIM];
}f_point_t;

//convert point_t to f_point_t
f_point_t f_point_t_convert(point_t point);

//convert point_t to f_point_t
point_t f_point_t_to_point(f_point_t p);

//check equility of two points
int f_point_t_eq(f_point_t p1, f_point_t p2);

typedef struct line_t{
	int pnt_cnt;
	point_t* line;
}line_t;

//debug print of line_t object
void line_t_dump(line_t line);

//return bary coordinates of triangle placed on "a", "b", "c"
point_t point_t_bary_coord( const point_t a,
                            const point_t b,
                            const point_t c,
                            const point_t p);

//return (p_2 - p_1)^\ortho
point_t get_ortho_vector(point_t p1, point_t p2, point_t p3);

#ifdef __cplusplus
}
#endif


#endif
