#ifndef _BOX_H
#define _BOX_H

#include "data.h"
#include "nets.h"

#ifdef __cplusplus
extern "C" {
#endif

#define DEFOLT_DATA_ELEM_SIZE 4

//coord(i, j, k, n) -> (i * n_j * n_k + j * n_k + k) * n_n + n
//coord(n, i, j, k) -> n * n_i * n_j * n_k + i * n_j * n_k + j * n_k + k

typedef struct box_t{
	point_t start;
	double step;
	int borders[3];
	int nets_cnt;
	int dyn_nets_cnt;
	int dyn_data_size;
	data_t data;
} box_t;

typedef struct crd_t{
	int crd[4];
} crd_t;
//data[coord][nets.count] - массив, где нулевое число - его alloc-размер (c учётом этого числа),
//а первое число - занятый размер (без учёта первых чисел)

//return coordinates of block which content point "coord" of net with id "net_id"
crd_t box_t_get_crd(box_t box, point_t coord, int net_id);

//print spatial coords of block which content point "coord"
void box_t_print_pnt_ids(box_t box, point_t pnt);

//check whether point belongs to box
int point_belong_box(point_t pnt, box_t box);

//return end point of box, point that opposite the start
point_t box_t_get_end(box_t box);

//check crd is belongs to box
int crd_is_acceptable(box_t box, crd_t crd);

//convert crd to coord
int crd_to_coord(box_t box, crd_t crd);

//init box correspondingly static nets
void init_static(nets_t nets, box_t* box);

//create box for "nets"
//step - set distance between neighbour box's blocks
box_t box_t_construct(nets_t nets, double step);
void box_t_destruct(box_t* box); //destructor

//extend box to contain point crd
void recreate_box(nets_t nets, box_t* box, crd_t crd);

//update box correspondingly coords of nets nodes
void nets_t_update_box(nets_t nets, box_t* box);

//sets exactly which box blocks the element's surface intersects with
//ATTENTION: expensive operation
void box_t_set_elem_exact(nets_t nets, box_t* box);

//debug print of crd_t object
void crd_t_dump(crd_t crd);

//debug print of box_t object
void box_t_dump(box_t box);

//debug print of box borders
void print_brds(double brds[3][2]);

//find projection of point to net if point is quite close to net
//if (sqr_distance != NULL) set square of distance from point to the projection here
void box_t_local_point_to_net_projection(box_t box, point_t point, net_t net, int net_id, double* sqr_distance);

//find projection of node coord to net if point is quite close to net
//if (sqr_distance != NULL) set square of distance from node to the projection here
void box_t_local_node_to_net_projection(box_t box, node_t* node, net_t net, int net_id, double* sqr_distance);

void extend_borders(double coef, double brds[3][2]);
void nets_t_fill_borders(nets_t nets, double res[3][2]);

#ifdef __cplusplus
}
#endif

#endif
