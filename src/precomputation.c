#include "precomputation.h"
#include <stdlib.h>

#ifndef M_PI
    #define M_PI 3.141592653589793238462643383279502884
#endif

double Contact_Resolution = 1.5;
double Contact_Const = 0.7;
double Min_Contact_Const = 0.5;

double Allow_shift = 0.5e-2 * 0.5e-2;
double Max_shift = 2e-2 * 2e-2;

point_t Fiber_Direction;
double** __Dihedral_Angles;
//############precomputation############################################
double max_3d(double x1, double x2, double x3){
	if (x1 > x2 && x1 > x3) return x1;
	if (x2 > x3) return x2;
	else return x3;
}

void auto_set_contact_recognition_consts(nets_t nets, double* min_resolution){
	double d = 0;
	unsigned int nets_cnt_elem = 0, cnt = nets.count, i;
	for (i = 0; i < cnt; i++){
		if (net_is_static(nets.nets[i])) continue;
		unsigned int n_elems = nets.nets[i].elems.count, j;
		for (j = 0; j < n_elems; j++){
			elem_t* elem = nets.nets[i].elems.elems[j];
			point_t cntr = (*elem).cntr_mass;
			double x1 = point_t_length(point_t_dif(cntr, (*((*elem).vrts[0])).coord));
			double x2 = point_t_length(point_t_dif(cntr, (*((*elem).vrts[1])).coord));
			double x3 = point_t_length(point_t_dif(cntr, (*((*elem).vrts[2])).coord));
			d += max_3d(x1, x2, x3);
		}
		nets_cnt_elem += n_elems;
	}
	d = d / nets_cnt_elem;
	if (min_resolution != NULL && d > *min_resolution)
		d = *min_resolution;
	double coef = 1.1;//1.05;
	Contact_Resolution = d * coef;
	//if (Contact_Resolution < Min_Contact_Const) Contact_Resolution = Min_Contact_Const; //experemental acceleration
	Contact_Resolution *= Contact_Resolution;
}

int net_t_count_shared_elems(node_t* node1, node_t* node2, int* elems_id){ //put into elems_id shared elems
	if (elems_id != NULL)
		elems_id[0] = -1, elems_id[1] = -1;
	unsigned int cnt1 = (*node1).cnt_elems, cnt2 = (*node2).cnt_elems;
	unsigned int i = 0, j = 0, k = 0;
	while (i < cnt1 && j < cnt2 && k < 2){
		if ((*node1).elems_id[i] < (*node2).elems_id[j]) ++i;
		else {
			if ((*node1).elems_id[i] > (*node2).elems_id[j]) ++j;
			else {
				if (elems_id != NULL) elems_id[k] = (*node1).elems_id[i];
				++k;
				++i;
			}
		}
	}
	return k;
}

void net_t_recognize_bnds(net_t net){
	unsigned int nv = net.vrtx.count, i;
	for (i = 0; i < nv; i++){
		node_t* origin = net.vrtx.nodes[i];
		int mob = 0, unmob = 0;
		unsigned int cnt_shar_spr = (*origin).cnt_springs, j;
		for (j = 0; j < cnt_shar_spr; j++){
			spring_t* spr = net.springs.springs[origin[0].springs_id[j]];
			node_t* neigh = spring_t_get_other_end(spr, origin);

			if (net_t_count_shared_elems(origin, neigh, NULL) < 2){
				if (is_free(neigh[0].state)) mob++;
				if (is_fix(neigh[0].state)) unmob++;
			}
		}
		if (is_fix(origin[0].state)){
			if (!mob) origin[0].state = FIX_BND;
			else if (mob && unmob) origin[0].state = EXTR_BND;
		}
		else {
			if (mob || unmob) origin[0].state = MOB_BND;
			else origin[0].state = IN;
		}
	}
}

void nets_t_recognize_bnds(nets_t nets){
	unsigned int i, cnt = nets.count;
	for (i = 0; i < cnt; i++)
		net_t_recognize_bnds(nets.nets[i]);
}

double spring_t_set_coef(net_t net, spring_t* obj){
	node_t* node1 = (*obj).ends[0];
	node_t*	node2 = (*obj).ends[1];
	int elems_id[2] = {};
	int cnt = net_t_count_shared_elems(node1, node2, elems_id);
	double area = net.elems.elems[elems_id[0]][0].area_0;
	if (cnt > 1)
		area += net.elems.elems[elems_id[1]][0].area_0;
	double h = ((*node1).h + (*node2).h) / 2;

	double coef = area * h / (*obj).l_0;
	return coef;
}


void spring_t_set_stress_params(spring_t* obj, double coef){
	point_t fiber_direction = Fiber_Direction;
	(*obj).stress_params = stress_t_construct((*obj).direction, &fiber_direction, coef);
	return;
}

void springs_t_set_params(net_t net){
	unsigned int cnt = net.springs.count, i;
	for (i = 0; i < cnt; i++){
		spring_t* spring = net.springs.springs[i];
		double coef = spring_t_set_coef(net, spring);
		spring_t_set_stress_params(spring, coef);
	}
}

int __TEMPVAR = 0;

void net_t_set_relax_state(net_t net){
	int e_cnt = net.elems.count, s_cnt = net.springs.count, i;
	for (i = 0; i < e_cnt; i++){
		elem_t* elem = net.elems.elems[i];
		elem[0].area_0 = point_t_length(elem[0].or_area);
	}
	for (i = 0; i < s_cnt; i++){
		spring_t* spr = net.springs.springs[i];
		spr[0].l_0 = spr[0].l;
	}
	springs_t_set_params(net);

	__Dihedral_Angles[__TEMPVAR] = (double*) calloc(s_cnt, sizeof(double));
	for (i = 0; i < s_cnt; ++i){
        spring_t* spr = net.springs.springs[i];
        node_t** ends = spr->ends;
        int se_id[2];
        int k = net_t_count_shared_elems(ends[0], ends[1], se_id);
        __Dihedral_Angles[__TEMPVAR][i] = 4 * M_PI;
        spr->isdigedral = 0;
        for (int l = 0; l < k; ++l)
        for (int j = 0; j < 3; ++j)
            if (net.elems.elems[se_id[l]]->vrts[j] != ends[0] && net.elems.elems[se_id[l]]->vrts[j] != ends[1]){
                spr->dihedral[l] = net.elems.elems[se_id[l]]->vrts[j];
                break;
            }
        if (k == 2) {
            spr->isdigedral = 1;
            point_t p[4];
            p[0] = spr->dihedral[0]->coord, p[1] = spr->dihedral[1]->coord;
            p[2] = ends[0]->coord, p[3] = ends[1]->coord;
            point_t n[2];
            for (int j = 0; j < 2; ++j)
                n[j] = NORM(point_t_or_area(p[2], p[3], p[j]));
            SCAL_S(-1, &n[1]);
            point_t res = CROSS(n[0], n[1]);
            double ddif = (1 - DOT(n[0], n[1])) / 2;
            if (ddif < 0) ddif = 0;
            if (ddif > 1) ddif = 1;
            __Dihedral_Angles[__TEMPVAR][i] = 2 * asin(sqrt(ddif));
            if (DOT(res, DIF(p[3], p[2])) < 0)
                __Dihedral_Angles[__TEMPVAR][i] *= -1;
        }
	}
}

void nets_t_set_relax_state(nets_t nets, point_t fiber_dir){
    __Dihedral_Angles = (double**) calloc(nets.count, sizeof(double*));
	Fiber_Direction = fiber_dir;
	unsigned int cnt = nets.count, j;
	for (j = 0; j < cnt; j++){
        __TEMPVAR = j;
		net_t_set_relax_state(nets.nets[j]);
    }
}

void net_t_elem_t_set_elems_neighbours(net_t net, elem_t* elem){
    node_t** vrtx = (*elem).vrts;
    unsigned int cnt1 = (*vrtx[0]).cnt_elems, cnt2 = (*vrtx[1]).cnt_elems, cnt3 = (*vrtx[2]).cnt_elems;
    int cnt_neighbours = cnt1 + cnt2 + cnt3;
    cnt_neighbours -= net_t_count_shared_elems(vrtx[0], vrtx[1], NULL);
    cnt_neighbours -= net_t_count_shared_elems(vrtx[0], vrtx[2], NULL);
    cnt_neighbours -= net_t_count_shared_elems(vrtx[1], vrtx[2], NULL);
    cnt_neighbours++;
    (*elem).cnt_neighbours = cnt_neighbours;
    (*elem).neighbours_id = (unsigned int*) calloc(cnt_neighbours, sizeof(unsigned int));
    unsigned int cur_id = (*elem).id;
    unsigned int i, k = 0;
    for (i = 0; i < cnt1; i++){
        unsigned int neigh_id = (*vrtx[0]).elems_id[i];
        if (neigh_id != cur_id)
            ((*elem).neighbours_id)[k++] = neigh_id;
    }
    for (i = 0; i < cnt2; i++){
        unsigned int neigh_id = (*vrtx[1]).elems_id[i];
        if (neigh_id != cur_id && !check_uint_in((*elem).neighbours_id, k, neigh_id))
            (*elem).neighbours_id[k++] = neigh_id;
    }
    for (i = 0; i < cnt3; i++){
        unsigned int neigh_id = (*vrtx[2]).elems_id[i];
        if (neigh_id != cur_id && !check_uint_in((*elem).neighbours_id, k, neigh_id))
            (*elem).neighbours_id[k++] = neigh_id;
    }
}

void net_t_set_elems_neighbours(net_t net){
	unsigned int e_cnt = net.elems.count, j;
	for (j = 0; j < e_cnt; j++){
		elem_t* elem = net.elems.elems[j];
        net_t_elem_t_set_elems_neighbours(net, elem);
	}
}

void nets_t_set_elems_neighbours(nets_t nets){
	unsigned int i, cnt = nets.count;
	for (i = 0; i < cnt; i++)
		net_t_set_elems_neighbours(nets.nets[i]);
}

void set_relax_consts(nets_t nets, double rel_allow_shift, double rel_max_shift){
	int nets_c = nets.count, i;
	double max = 0, min = nets.nets[0].springs.springs[0][0].l_0;
	for(i = 0; i < nets_c; i++){
		net_t net = nets.nets[i];
		if (net_is_static(net)) continue;
		int s_c = net.springs.count, j;
		for (j = 0; j < s_c; j++){
			double len = net.springs.springs[j][0].l_0;
			if (max < len) max = len;
			if (min > len) min = len;
		}
	}
	Allow_shift = min * rel_allow_shift;
	Allow_shift *= Allow_shift;
	Max_shift = min * rel_max_shift;
	Max_shift *= Max_shift;
}

void precomputation(nets_t nets){
	//init_nets(nets);
	nets_t_recognize_bnds(nets);
	//nets_t_set_springs_params(nets);		//order is important
	nets_t_construct_nodes_contact(nets);	//order is unimportant
	nets_t_set_elems_neighbours(nets);		//order is unimportant
	set_relax_consts(nets, 0.02, 0.05);		//0.02, 0.05
	auto_set_contact_recognition_consts(nets, NULL);	//order is important
}
//############precomputation############################################

void set_contact_recognition_resolution(double contact_resolution){
	Contact_Resolution = contact_resolution * contact_resolution;
}

void set_contact_to_plate_recognition_const(double contact_const){
	Contact_Const = contact_const * contact_const;
}

void set_min_contact_recognition_const(double min_contact_const){
	Min_Contact_Const = min_contact_const;
}

void set_contact_recognition_consts(double contact_resolution, double contact_const, double min_contact_const){
	set_contact_recognition_resolution(contact_resolution);
	set_contact_to_plate_recognition_const(contact_const);
	set_min_contact_recognition_const(min_contact_const);
}

double get_Contact_Resolution(){
	return sqrt(Contact_Resolution);
}

double get_Contact_Const(){
	return sqrt(Contact_Const);
}
double get_Min_Contact_Const(){
	return Min_Contact_Const;
}
