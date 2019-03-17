#ifndef _INTERSECT_H
#define _INTERSECT_H

#include "nets.h"
#include "bound-box.h"

//return 0 if there is no intersection between line and "bar" net of nets
//box - aligned bounded box around nets, it must be initialized with "box_t_set_elem_exact" for proper work of function
int line_to_boxed_nets_intersection(point_t line[2], nets_t nets, int bar, box_t box, point_t* intersect);

#endif 
