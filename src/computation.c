#include <float.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include "save-data.h"
#include "format_in.h"
#include "precomputation.h"
#include "computation.h"

#ifndef M_PI
    #define M_PI 3.141592653589793238462643383279502884
#endif

double Contact_Force_Const1 = 0.0005;//0.001;//0.00008;
double Contact_Force_Const2 = 0.02;//0.01;//0.02;
int ElasticModelType = 0;

void set_contact_force_consts(double k1, double k2){
	Contact_Force_Const1 = k1;
	Contact_Force_Const1 = k2;
}
//############computation###############################################
point_t pressure_force(net_t net, node_t* node, double P){
	unsigned int e_cnt = (*node).cnt_elems, i;
	point_t f = get_zero_point();
	for (i = 0; i < e_cnt; i++){
		f = point_t_sum(f, (*(net.elems.elems[(*node).elems_id[i]])).or_area);
	}
	point_t_coef_mul(P/3 * 133.3, &f); //10^-6 * H
	return f;
}

double get_min_len(net_t net, node_t* node){
	int s_cnt = node[0].cnt_springs, i;
	double min_len = net.springs.springs[node[0].springs_id[0]][0].l_0;
	for (i = 0; i < s_cnt; i++){
		double l = net.springs.springs[node[0].springs_id[i]][0].l_0;
		if (min_len > l) min_len = l;
	}

	return min_len;
}

double recommended_delta(nets_t nets, double P){//works not very good
	int nets_cnt = nets.count, j, k = 0;
	double eps = 0.005;
	double res1 = 20;//0;//20;
	double res2 = 0;
	for (j = 0; j < nets_cnt; j++){
		net_t net = nets.nets[j];
		if (net_is_static(net)) continue;
		int n_cnt = net.vrtx.count, i;
		for (i = 0; i < n_cnt; i++){
			node_t* node = net.vrtx.nodes[i];
			point_t F_p = pressure_force(net, node, P);
			double f_p = point_t_length(F_p);
			double l = get_min_len(net, node);
			double delta = eps * l / f_p;
			if (res1 > delta) res1 = delta;
			res2 += delta;
			k++;
		}
	}
	res2 /= k;
	return (res1 + res2) / 2;
}

point_t f_spring(net_t net, node_t* node, spring_t* spring){
	point_t force = spring[0].direction;
	SCAL_S(((spring[0].ends[1] == node) ? -1 : 1) * spring[0].force, &force);

	return force;
}


point_t elastic_force_MSM(net_t net, node_t* node){
	unsigned int sp_cnt = (*node).cnt_springs, i;
	point_t f = get_zero_point();
	for (i = 0; i < sp_cnt; i++){
		f = point_t_sum(f, f_spring(net, node, net.springs.springs[(*node).springs_id[i]]));
	}

	return f; //10^-4 * H
}

point_t elastic_force_REINFORCING(net_t net, node_t* node){
    double mu = 1000000.0  / 3.0;
    double Jm = 2.3;

    point_t f = ZERO();

    for (unsigned int i = 0; i < node->cnt_elems; ++i)
    {
        elem_t* triag = net.elems.elems[node->elems_id[i]];
        node_t* n[3] = {};
        int flag = 0;
        for (int j = 0; j < 3; ++j)
        {
            n[j] = triag->vrts[j];
            if (triag->vrts[j] == node)
            {
                n[j] = n[0];
                n[0] = node;
                flag = 1;
            }
        }
        assert(flag);

        point_t ff = ZERO();
        assert(LEN(ff) < 10 * DBL_EPSILON);
        double Ap = triag->area_0;
        point_t D[3] = {};
        for (int j = 0; j < 3; ++j)
            D[j] = SCAL(1.0 / 2 / Ap, get_ortho_vector(n[(j + 1)%3]->initial, n[(j + 2)%3]->initial, n[(j)%3]->initial));
        double Aq = LEN(triag->or_area);

        double trC = 0;
        for (int j = 0; j < 3; ++j)
            for (int k = 0; k < j; ++k)
                trC += 2 * DOT(n[j]->coord, n[k]->coord) * DOT(D[j], D[k]);
        for (int j = 0; j < 3; ++j)
            trC += DOT(n[j]->coord, n[j]->coord) * DOT(D[j], D[j]);

        double J = Aq / Ap;
        for (int j = 0; j < 3; ++j)
            ADD_S(&ff, Ap * DOT(D[0], D[j]), n[j]->coord);
        ADD_S(&ff, -1.0 / 2 / J / J / J, get_ortho_vector(n[1]->coord, n[2]->coord, n[0]->coord));
        ADD_S(&f, mu / (1 - (trC + 1 / J / J - 3) / Jm), ff);
    }

    SCAL_S(-node->h, &f);

    return f;
}

point_t elastic_force_NEOGOOK(net_t net, node_t* node){
    double mu = 1000000.0 / 3;
    //double d = 1;

    point_t f = ZERO();
    for (unsigned int i = 0; i < node->cnt_elems; ++i)
    {
        elem_t* triag = net.elems.elems[node->elems_id[i]];
        node_t* n[3] = {};
        int flag = 0;
        for (int j = 0; j < 3; ++j)
        {
            n[j] = triag->vrts[j];
            if (triag->vrts[j] == node)
            {
                n[j] = n[0];
                n[0] = node;
                flag = 1;
            }
        }
        assert(flag);

        double Ap = triag->area_0;
        point_t D[3] = {};
        for (int j = 0; j < 3; ++j)
            D[j] = SCAL(1.0 / 2.0 / Ap, get_ortho_vector(n[(j + 1)%3]->initial, n[(j + 2)%3]->initial, n[(j)%3]->initial));
        double Aq = LEN(triag->or_area);
//        double trC = 0;
//        for (int j = 0; j < 3; ++j)
//            for (int k = 0; k < j; ++k)
//                trC += 2 * DOT(n[j]->coord, n[k]->coord) * DOT(D[j], D[k]);
//        for (int j = 0; j < 3; ++j)
//            trC += DOT(n[j]->coord, n[j]->coord) * DOT(D[j], D[j]);
//        static int tmp = 0;
//        tmp++;
//        if (!(tmp % 100000)) printf("trC = %lg\n", trC);

        double J = Aq / Ap;
       // double Jm = 2.3;
       // double coef = 1.0 / (1 - (trC + 1 / J / J - 3) / Jm);
        for (int j = 0; j < 3; ++j)
            ADD_S(&f, mu * Ap * DOT(D[0], D[j]), n[j]->coord);

        ADD_S(&f, -mu / 2 / J / J / J, get_ortho_vector(n[1]->coord, n[2]->coord, n[0]->coord));
        //ADD_S(&f, ((J - 1) * d - 1)* mu / 2, get_ortho_vector(n[1]->coord, n[2]->coord, n[0]->coord));
        /*for (int j = 0; j < 3; ++j)
            ADD_S(&f, mu * Ap / 2 * DOT(D[0], D[j]) * (trC - 3), n[j]->coord);
        ADD_S(&f, (J - 1) * d, get_ortho_vector(n[1]->initial, n[2]->initial, n[0]->initial));*/
    }

    SCAL_S(-node->h, &f);

    return f;

}



point_t elastic_force_TRQS(net_t net, node_t* node){     //St Venan
    double E = 1000000;      //Young modulus, Pa
	double nu = 0.45;           //Poisson ratio

    point_t f = ZERO();
    for (unsigned int i = 0; i < node->cnt_elems; ++i)
    {
        elem_t* triag = net.elems.elems[node->elems_id[i]];
        node_t* n[3] = {};
        int flag = 0;
        for (int j = 0; j < 3; ++j)
        {
            n[j] = triag->vrts[j];
            if (triag->vrts[j] == node)
            {
                n[j] = n[0];
                n[0] = node;
                flag = 1;
            }
        }
        assert(flag);

        point_t dir[3] = {};
        dir[0] = DIF(n[1]->coord, n[2]->coord);
        dir[1] = DIF(n[2]->coord, n[0]->coord);
        dir[2] = DIF(n[0]->coord, n[1]->coord);
        double dir_len[3] = {};
        for (int j = 0; j < 3; ++j)
            dir_len[j] = LEN(dir[j]);

        double k[3] = {}, c[3] = {}, alpha[3] = {};
        for (int j = 0; j < 3; ++j)
            alpha[j] = M_PI - acos (DOT(dir[(j + 2) % 3], dir[(j + 1) % 3]) / (dir_len[(j + 2) % 3] * dir_len[(j + 1) % 3]));
//        static int counter = 0;
//        if (counter == 16)
//            printf("smth");
//        if (n[0]->id == 200) counter++;
//        if (!(fabs(M_PI - alpha[0] - alpha[1] - alpha[2]) < 1.0e-6))
//        {
//
//            printf ("dif = %lg\n", fabs(M_PI - alpha[0] - alpha[1] - alpha[2]));
//            point_t_dump(n[0]->coord);
//            point_t_dump(n[0]->next);
//            printf("id = %d, cnt = %d\n", n[0]->id, counter++);
//            point_t_dump(n[1]->coord);
//            point_t_dump(n[2]->coord);
//        }
//        assert(fabs(M_PI - alpha[0] - alpha[1] - alpha[2]) < 1.0e-6);

        double ctg[3] = {};
        for (int j = 0; j < 3; ++j)
            ctg[j] = 1.0 / tan(alpha[j]);
        double scl = 16 * triag->area_0 * (1 - nu * nu) / E;
        for (int j = 0; j < 3; ++j)
        {
            k[j] = (2 * ctg[j] * ctg[j] + 1 - nu) / scl;
            c[j] = (2 * ctg[(j + 1) % 3] * ctg[(j + 2) % 3] + nu - 1) / scl;
        }

        double d2l[3] = {};
        for (int j = 0; j < 3; ++j)
        {
            spring_t* spr = get_shared_spring(net, n[(j + 1) % 3], n[(j + 2) % 3]);
            assert(spr != NULL);
            d2l[j] = dir_len[j] * dir_len[j] - spr->l_0 * spr->l_0;
        }
        SCAL_S(-1, &dir[2]);

        for (int j = 1; j < 3; ++j)
        {
            double cf = k[j] * d2l[j] + c[3 - j] * d2l[0] + c[0] * d2l[3 - j];

            ADD_S(&f, cf, dir[j]);
        }
    }

    SCAL_S(node->h, &f);
    //if (LEN(f) > 10000){
    //    printf("Foo");
    //}

    return f;   //10^-6 * H
}


/*
point_t elastic_force1(net_t net, node_t* node){
    double E = 1000000;         //Young modulus, Pa

    point_t f = ZERO();
    for (unsigned int i = 0; i < node->cnt_elems; ++i)
    {
        elem_t* triag = net.elems.elems[node->elems_id[i]];
        node_t* n[3];
        for (int j = 0; j < 3; ++j)
        {
            n[j] = triag->vrts[j];
            if (triag->vrts[j] == node)
            {
                n[j] = n[0];
                n[0] = node;
            }
        }

        point_t dir[3];
        dir[0] = DIF(n[1]->coord, n[2]->coord);
        dir[1] = DIF(n[2]->coord, n[0]->coord);
        dir[2] = DIF(n[0]->coord, n[1]->coord);
        double dir_len[3];
        for (int j = 0; j < 3; ++j)
            dir_len[j] = LEN(dir[j]);


        double kk = E * node->h;
        double area[3] = {};

        double dl[3] = {};
        double l_0[3] = {};
        spring_t* sprr[3];
        for (int j = 0; j < 3; ++j)
        {
            spring_t* spr = get_shared_spring(net, n[(j + 1) % 3], n[(j + 2) % 3]);
            assert(spr != NULL);
            dl[j] = dir_len[j]  -  spr->l_0;
            if (fabs(dir_len[j] - spr->l) > 1.0e-8)
                printf("str 195\n");
            l_0[j] = spr->l_0;
            double lll = spr->l;
            int elems_id[2];
            for (int k = 0, cnt = net_t_count_shared_elems(n[(j + 1) % 3], n[(j + 2) % 3], elems_id); k < cnt; ++k)
                area[j] += net.elems.elems[elems_id[k]]->area_0;
            sprr[j] = spr;
        }
        for (int j = 0; j < 3; ++j){
        double coef = sprr[j]->stress_params.sigma1 - E * node->h / l_0[j] * area[j];
        if (coef > 1e-8)
            printf("Странно\n");
        }

        SCAL_S(-1, &(dir[2]));

        for (int j = 1; j < 3; ++j)
        {
            double cf = kk * triag->area_0 / l_0[j] / l_0[j] * dl[j];
            NORM_S(&dir[j]);
            ADD_S(&f, cf, dir[j]);
            res[2 * i + j - 1] = SCAL(cf, dir[j]);
        }
    }
    for (int i = 0; i < node->cnt_springs; ++i)
    {
        double E = 1000000;
        spring_t* spr = net.springs.springs[node->springs_id[i]];
        int elems_id[2];
        double area = 0;
        for (int k = 0, cnt = net_t_count_shared_elems(node, spring_t_get_other_end(spr, node), elems_id); k < cnt; ++k)
                area += net.elems.elems[elems_id[k]]->area_0;
        //pes[3 + i%3] = SCAL(E * node->h / spr->l_0 * area * (spr->l - spr->l_0) / spr->l_0, spr->direction);
        node_t* node1 = spring_t_get_other_end(spr, node);
        point_t dir = DIF(node1->coord, node->coord);
        point_t dir1 = SCAL(spr->l, spr->direction);
        double dif = LEN(DIF(dir, dir1));
            if (LEN(SUM(dir, dir1)) < dif) dif = LEN(SUM(dir, dir1));
        if (dif > 1.0e-5)
            printf("Вот оно\n");
        double hhhh = 4;
    }


    return f;   //10^-6 * H
}*/

void set_elastic_model(int type) { ElasticModelType = type; }

point_t elastic_force(net_t net, node_t* node){
    typedef point_t (*elast_t)(net_t net, node_t* node);
    static elast_t elast_model[] = {
        elastic_force_MSM,
        elastic_force_TRQS,
        elastic_force_NEOGOOK,
        elastic_force_REINFORCING
    };

    return elast_model[ElasticModelType](net, node);
}

void get_foces_to_bending_elem(spring_t* spr, point_t f[4], double teta_0){
    double k_e = 100;

    node_t** ends = spr->ends;
    node_t** dihed = spr->dihedral;
    point_t n[2];
    double s[2];
    for (int j = 0; j < 2; ++j){
        n[j] = point_t_or_area(ends[0]->coord, ends[1]->coord, dihed[j]->coord);
        s[j] = LEN(n[j]);
        NORM_S(&n[j]);
    }

    SCAL_S(-1, &n[1]);
    point_t res = CROSS(n[0], n[1]);
    double ddif = (1 - DOT(n[0], n[1])) / 2;
    if (ddif < 0) ddif = 0;
    if (ddif > 1) ddif = 1;
    double teta = 2 * asin(sqrt(ddif));
    if (DOT(res, DIF(ends[1]->coord, ends[0]->coord)) < 0)
        teta *= -1;
    point_t p[4] = {};
    p[0] = spr->dihedral[0]->coord, p[1] = spr->dihedral[1]->coord;
    p[2] = ends[0]->coord, p[3] = ends[1]->coord;

    point_t E = DIF(p[3], p[2]);
    double e = LEN(E);
    NORM_S(&E);

    f[0] = SCAL(e / s[0], n[0]);
    f[1] = SCAL(e / s[1], n[1]);
    f[2] = ADD(SCAL(DOT(DIF(p[0], p[3]), E) / s[0], n[0]), DOT(DIF(p[1], p[3]), E) / s[1], n[1]);
    f[3] = ADD(SCAL(-DOT(DIF(p[0], p[2]), E) / s[0], n[0]), -DOT(DIF(p[1], p[2]), E) / s[1], n[1]);

    double coef = k_e * e * e / (s[0] + s[1]);
    for (int i = 0; i < 4; ++i)
        SCAL_S(coef * (sin(teta / 2) - sin(teta_0 / 2)), &f[i]); //sin(teta / 2) - sin(teta_0 / 2)
}

point_t* __bendings;
int __sz_bendings = 0;
int __TEMPVNET = 0;

void init_bending_forces(net_t net){
    if (__sz_bendings == 0){
        __bendings = (point_t*) calloc(net.vrtx.count, sizeof(point_t));
        __sz_bendings = net.vrtx.count;
    }
    else if (__sz_bendings < net.vrtx.count) {
        __bendings = (point_t*) realloc(__bendings, net.vrtx.count * sizeof(point_t));
        __sz_bendings = net.vrtx.count;
    }
    memset(__bendings, 0, __sz_bendings * sizeof(point_t));

    for (int i = 0, s_cnt = net.springs.count; i < s_cnt; ++i)
    {
        if (!net.springs.springs[i]->isdigedral) continue;
        point_t p[4] = {};
        spring_t* spr = net.springs.springs[i];
        get_foces_to_bending_elem(spr, p, __Dihedral_Angles[__TEMPVNET][i]);
        ADD_S(&__bendings[spr->dihedral[0]->id], 1, p[0]);
        ADD_S(&__bendings[spr->dihedral[1]->id], 1, p[1]);
        ADD_S(&__bendings[spr->ends[0]->id], 1, p[2]);
        ADD_S(&__bendings[spr->ends[1]->id], 1, p[3]);
    }

}

point_t bending_force(net_t net, node_t* node, int vrt_id)
{
    return __bendings[vrt_id];;
}

point_t get_contact_force(double dist, point_t normal, double force, node_t* node){
	double k1 = Contact_Force_Const1;
	double k2 = Contact_Force_Const2;
	double coef = 0;
	if (dist > 0)
		coef = force * exp((-1) * k1 * dist / force);
	else
		coef = force + k2 * (-dist);
	point_t f = SCAL(coef, normal);
	return f;
}

int check_contact(elem_t* elem, node_t* node, point_t* f, double* r){
	point_t cur_shift = point_t_dif((*node).coord, (*elem).cntr_mass);
	double node_elem_cntr_dist = DOT(cur_shift, cur_shift);
	if (node_elem_cntr_dist > Contact_Resolution) return 0;

	point_t normal = NORM((*elem).or_area);
	double plane_dist = DOT(DIF((*node).next, (*elem).cntr_mass), normal);
	if (plane_dist < Contact_Const) {
		(*node).coaptative = 1;
		if (r != NULL) *r = node_elem_cntr_dist;
		double force = LEN(DIF((*node).next, (*node).coord));
		*f = SUM(*f, get_contact_force(plane_dist, normal, force, node));

		return 1;
	}
	return 0;
}

int check_locale_contact_neigh_elem(nets_t nets, int bar, node_t* node, point_t* f, int set_new_id, double r){
	if ((*node).contact_elem_id == NULL || (*node).contact_elem_id[bar] < 0) return 0;
	elem_t* loc_elem = nets.nets[bar].elems.elems[(*node).contact_elem_id[bar]];
	int ret = 0;
	int i, cnt_neighbours = (*loc_elem).cnt_neighbours;
	if (!set_new_id) r = 1000;
	for (i = 0; i < cnt_neighbours; i++){
		unsigned int neigh_id = (*loc_elem).neighbours_id[i];
		elem_t* neighbour = nets.nets[bar].elems.elems[neigh_id];
		double cur_r = 0;
		int curret = check_contact(neighbour, node, f, &cur_r);
		if ((!set_new_id || cur_r < r) && curret ){
			(*node).contact_elem_id[bar] = neigh_id;
			r = cur_r;
		}
		ret += curret;
	}
	if (!ret && set_new_id) (*node).contact_elem_id[bar] = -1;
	return ret;
}

int compute_locale_contact_elems(nets_t nets, int bar, node_t* node, point_t* f){
	if ((*node).contact_elem_id == NULL || (*node).contact_elem_id[bar] < 0) return 0;
	elem_t* loc_elem = nets.nets[bar].elems.elems[(*node).contact_elem_id[bar]];
	double r = 0;
	int ret = check_contact(loc_elem, node, f, &r);
	ret += check_locale_contact_neigh_elem(nets, bar, node, f, !ret, r);
	if (!ret) (*node).contact_elem_id[bar] = -1;
	else point_t_coef_mul(1.0 / ret, f);
	return ret;
}

point_t contact_force(nets_t nets, int orig, int bar, node_t* node){
	net_t barrier = nets.nets[bar];
	unsigned int e_cnt = barrier.elems.count, i;
	point_t f = get_zero_point();
	if (compute_locale_contact_elems(nets, bar, node, &f)) return f;

	int cnt = 0;
	double r = 1000;
	for (i = 0; i < e_cnt; i++){
		point_t cur_f = get_zero_point();
		elem_t* elem = barrier.elems.elems[i];
		double cur_r = 0;
		if (check_contact(elem, node, &cur_f, &cur_r)) {
			if (cur_r < r){
				(*node).contact_elem_id[bar] = i;
				r = cur_r;
				f = cur_f;
			}
			cnt = 1;
		}
	}
	//if (cnt) printf("orig_net: %d node: %d bar_net: %d bar_node: %d\n", orig, node[0].id, bar, (*node).contact_elem_id[bar]);
	cnt += check_locale_contact_neigh_elem(nets, bar, node, &f, 0, 0);
	if (cnt) point_t_coef_mul(1.0 / cnt, &f);
	return f;
}


#define CHECK_ACCEPTABLE(X, ID)\
curcrd.crd[X] = crd.crd[X] + ID;\
if (!(curcrd.crd[X] >= 0 && curcrd.crd[X] < box.borders[X]))\
	continue

point_t box_contact_force(nets_t nets, int orig, int bar, node_t* node, box_t box){
	point_t f = get_zero_point();
	if (compute_locale_contact_elems(nets, bar, node, &f)) return f;

	net_t barrier = nets.nets[bar];
	data_t data = box.data;
	crd_t crd = box_t_get_crd(box, node[0].coord, bar); //it's possible to optimize
	//printf("gr");
	//crd_t_dump(crd);
	int cnt = 0;
	double r = 1000;
	int i, j, k;
	crd_t curcrd;
	curcrd.crd[3] = bar;
	for (i = -1; i <= 1; i++){ CHECK_ACCEPTABLE(0, i);
	for (j = -1; j <= 1; j++){ CHECK_ACCEPTABLE(1, j);
	for (k = -1; k <= 1; k++){ CHECK_ACCEPTABLE(2, k);
		int coord = crd_to_coord(box, curcrd);
		unsigned int e_cnt = data.busy_len[coord];
		unsigned int l;
		for (l = 0; l < e_cnt; l++){
			int elem_id = data.data[coord][l];
			elem_t* elem = barrier.elems.elems[elem_id];
			point_t cur_f = get_zero_point();
			double cur_r = 0;
			if (check_contact(elem, node, &cur_f, &cur_r)) {
			//if (check_contact(elem, node, &f, &cur_r)) {
				if (cur_r < r){
					(*node).contact_elem_id[bar] = elem_id;
					r = cur_r;
					f = cur_f;
				}
				cnt = 1;
			}
		}
	}}}
	//if (cnt) printf("orig_net: %d node: %d bar_net: %d bar_node: %d\n", orig, node[0].id, bar, (*node).contact_elem_id[bar]);
	cnt += check_locale_contact_neigh_elem(nets, bar, node, &f, 0, 0);
	if (cnt) point_t_coef_mul(1.0 / cnt, &f);
	return f;
}
#undef CHECK_ACCEPTABLE

void relaxation(nets_t nets, double coef){
	int cnt = nets.count, j;
	for (j = 0; j < cnt; j++){
		net_t net = nets.nets[j];
		if (net_is_static(net)) continue;
		int cnt_vrt = net.vrtx.count, i;
		for (i = 0; i < cnt_vrt; i++){
			node_t* node = net.vrtx.nodes[i];
			point_t shift = SCAL(coef, DIF((*node).next, (*node).coord));
			point_t next = SUM((*node).coord, shift);
			point_t_cpy_points(&next, &(*node).next);
		}
	}
}

double gt_elastic = 0;

double compute_free_nexts(net_t net, double P, double delta){
	int cnt_vrt = net.vrtx.count, i = 0;
	double max_shift = 0;
	//init_bending_forces(net);                     //bending
	for (i = 0; i < cnt_vrt; i++){
		node_t* node = net.vrtx.nodes[i];
		if (is_fix((*node).state)) continue;
		point_t F_p = pressure_force(net, node, P);

		struct timeval start, end;
        gettimeofday(&start, NULL);
		point_t F_sp = elastic_force(net, node);
		gettimeofday(&end, NULL);
        gt_elastic += (end.tv_sec  - start.tv_sec) * 1.e+3 + (end.tv_usec - start.tv_usec) * 1.e-3;

		//point_t F_bn = bending_force(net, node, i);       //bending
		//printf("p = "); point_t_dump(F_p);
		//printf("bn = "); point_t_dump(F_bn);
		//printf("sp = "); point_t_dump(F_sp);
		//ADD_S(&F_sp, 1.0, F_bn);

		point_t next = SCAL(delta, SUM(F_p, F_sp));

		double next_sqr_len = SQR_LEN(next);
		if (max_shift < next_sqr_len) max_shift = next_sqr_len;

		next = SUM(next, (*node).coord);
		point_t_cpy_points(&next, &(*node).next);
		(*(net.vrtx.nodes[i])).coaptative = -1;
	}
	return max_shift;
}

double get_diviation(nets_t nets){
	double diviation = 0;
	unsigned int nets_cnt = nets.count, i, j;
	for (i = 0; i < nets_cnt; i++){
		if (net_is_static(nets.nets[i])) continue;
		unsigned int node_cnt = nets.nets[i].vrtx.count;
		for (j = 0; j < node_cnt; j++){
			node_t* node = nets.nets[i].vrtx.nodes[j];
			point_t shift = DIF((*node).next, (*node).coord);
			diviation += SQR_LEN(shift);
		}
	}
	diviation = sqrt(diviation);
	return diviation;
}

double compute_constraint_nexts(nets_t nets, int flag, box_t box){
	unsigned int nets_cnt = nets.count;
	if (nets_cnt < 2) return -1;

	for (unsigned int i = 0; i < nets_cnt; i++){
		net_t origin = nets.nets[i];
		if (net_is_static(origin)) continue;
		unsigned int node_cnt = origin.vrtx.count;
		for (unsigned int k = 0; k < node_cnt; k++){
			node_t* node = origin.vrtx.nodes[k];
			if (is_fix((*node).state)) continue;
			point_t f_contact = get_zero_point();
			int stat = 0;
			unsigned int j = 0;
			for (j = 0; j < nets_cnt; j++){
				if (j == i) continue;
				if (!flag || (*node).contact_elem_id[j] >= 0){
					//point_t f_contact = contact_force(nets, i, j, node);
					f_contact = point_t_sum(f_contact, box_contact_force(nets, i, j, node, box));
					stat += ((*node).contact_elem_id[j] >= 0);
				}
			}
			if (stat){
				point_t next = point_t_sum((*node).next, f_contact);
				point_t_cpy_points(&next, &(*node).next);
			}
		}
	}
	return 0;//get_diviation(nets);
}

int compute_nexts(nets_t nets, double P, double delta, int constraint, box_t box){
	unsigned int nets_cnt = nets.count, i;
	double shift = 0;
	int ret = 0;
	for (i = 0; i < nets_cnt; i++){
		if (net_is_static(nets.nets[i])) continue;
		__TEMPVNET = i;
		double max_net_shift = compute_free_nexts(nets.nets[i], P, delta);
		if (shift < max_net_shift) shift = max_net_shift;
	}

	if (shift > Max_shift && FLAG_REL) {
		relaxation(nets, sqrt(Allow_shift / shift));
		ret++;
	}
	compute_constraint_nexts(nets, constraint, box);
	return ret;
}

point_t* Mid_Shift;
unsigned int Cnt_Nodes;

unsigned int get_Cnt_Nodes() { return Cnt_Nodes; }
point_t get_Mid_Shift(int i) { return Mid_Shift[i]; }
void construct_Mid_Shift(nets_t nets){
	unsigned int cnt_nodes = 0, cnt_nets = nets.count, i;
	for (i = 0; i < cnt_nets; i++)
		if (!net_is_static(nets.nets[i])) cnt_nodes += nets.nets[i].vrtx.count;
	Cnt_Nodes = cnt_nodes;
	Mid_Shift = (point_t*) calloc(cnt_nodes, sizeof(point_t));
}

void restart_Mid_Shift(){
	unsigned int i;
	for (i = 0; i < Cnt_Nodes; i++)
		Mid_Shift[i] = get_zero_point();
}

void update_Mid_Shift(nets_t nets){
	unsigned int nets_cnt = nets.count, i, j, cnt = 0;
	for (i = 0; i < nets_cnt; i++){
		if (net_is_static(nets.nets[i])) continue;
		unsigned int node_cnt = nets.nets[i].vrtx.count;
		for (j = 0; j < node_cnt; j++){
			node_t* node = nets.nets[i].vrtx.nodes[j];
			point_t shift = point_t_dif((*node).next, (*node).coord);
			Mid_Shift[cnt] = point_t_sum(Mid_Shift[cnt], shift);
			cnt++;
		}
	}
}

double get_Mid_Shift_len(int nsteps){
	unsigned int i;
	double res = 0;
	for (i = 0; i < Cnt_Nodes; i++)
		res += SQR_LEN(Mid_Shift[i]);
	res = sqrt(res);
	res /= nsteps;
	return res;
}

void compute_nets(nets_t nets, double P, double delta, unsigned int n, double eps, int freq){

	box_t box = box_t_construct(nets, get_Contact_Resolution());
	//box_t_dump(box);

	construct_Mid_Shift(nets);
	unsigned int i, k = 0;
	double init_div = 0;
	double res = 0;
	double constr_freq = (sqrt(Contact_Resolution / Max_shift) - 1);
	int step = (constr_freq > 10) ? (int) constr_freq : 10;
	printf("constr step = %d\n", step);
	int crush = 0;
	for (i = 1; i <= n && k < 3; i++){
		nets_t_update_box(nets, &box);

		int flag = i % 1000, constraint = i % step;

		crush += compute_nexts(nets, P, delta, constraint, box);

		if (!flag){
			 double diviation = get_diviation(nets);
			 printf("it = %d, diviation = %lg, relation = %lg\n", i, diviation / delta, diviation / init_div);
		 }

		update_Mid_Shift(nets);

		if (i == 1) {
			double max_Mid_Shift = 0;
			unsigned int l;
			for (l = 0; l < Cnt_Nodes; l++){
				double len = SQR_LEN(Mid_Shift[l]);
				if (max_Mid_Shift < len) max_Mid_Shift = len;
			}
			printf("Max_Shift = %e\n", sqrt(max_Mid_Shift));
			init_div = get_diviation(nets);
		}
		if (!(i % freq)){
			res = get_Mid_Shift_len(freq);
			double rms = res / Cnt_Nodes; // delta ;
			printf("RMS shift per iter of node = %e mm, relation = %lg, cr = %d\n", rms , res / init_div, crush);
			crush = 0;
			restart_Mid_Shift();
			if (res < eps * init_div) k++;
				else k = 0;
		}

		update_nets(nets);
	}
	printf("Made %u iterations, relation = %lg\n", --i, res / init_div);

	box_t_destruct(&box);
	free(Mid_Shift);
}



world_t* world_t_construct(nets_t dynamic_nets, nets_t static_nets, solver_t solver_data, collision_t collision_data, wrld_cnd_t conditions){
    world_t* world = (world_t*)calloc(1, sizeof(world_t));
    world->dynamic_nets = dynamic_nets;
    world->static_nets = static_nets;
    world->collision = collision_data;
    world->conditions = conditions;
    world->solver_data = solver_data;
    world->statistic = statistical_data_construct(dynamic_nets);
    world->union_nets = create_union_net(dynamic_nets, static_nets);
    set_initial_solving_params(world);
    return world;
}

void world_t_destruct(world_t* world){
    statistical_data_destruct(&world->statistic);
    nets_t_surfacial_free(world->union_nets);
    free(world);
}

statistic_t statistical_data_construct(nets_t nets){
    unsigned int cnt_nodes = 0, cnt_nets = nets.count, i;
	for (i = 0; i < cnt_nets; i++)
        cnt_nodes += nets.nets[i].vrtx.count;
	point_t* mid_shift = (point_t*) calloc(cnt_nodes, sizeof(point_t));
    statistic_t stat = {mid_shift, cnt_nodes, 0, 1.0};
    return stat;
}

void statistical_data_destruct(statistic_t *st){
    st->cnt_nodes = 0;
    free(st->mid_shift);
}

nets_t create_union_net(nets_t nets1, nets_t nets2){
    unsigned int un_cnt = nets1.count + nets2.count;
    nets_t nets = nets_t_get_net(un_cnt);
    for (unsigned int i = 0; i < nets1.count; ++i)
        nets.nets[i] = nets1.nets[i];
    for (unsigned int i = 0; i < nets2.count; ++i)
        nets.nets[i + nets1.count] = nets2.nets[i];

    return nets;
}

collision_t collision_data_t_construct(nets_t dynamic_nets, nets_t static_nets, int check_freq){
    double constr_freq = (sqrt(Contact_Resolution / Max_shift) - 1);
	int freq = (constr_freq > check_freq) ? (int) constr_freq : check_freq;
    nets_t union_nets = create_union_net(dynamic_nets, static_nets);
	box_t box = box_t_construct(union_nets, get_Contact_Resolution());
	nets_t_surfacial_free(union_nets);
	collision_t collision = {box, freq, -1};
	return collision;
}

void collision_data_t_destruct(collision_t cl){
    box_t_destruct(&cl.box);
    cl.check_freq = -1;
    cl.counter = -1;
}

void collision_data_t_update(nets_t nets, collision_t* coll_data){
    nets_t_update_box(nets, &coll_data->box);
    ++coll_data->counter;
}

int collision_t_get_state(collision_t coll_data){
    return coll_data.counter % coll_data.check_freq;
}

void update_statistic_t(nets_t nets, statistic_t *st){
    unsigned int nets_cnt = nets.count, cnt = 0;
	for (unsigned int i = 0; i < nets_cnt; ++i){
		unsigned int node_cnt = nets.nets[i].vrtx.count;
		for (unsigned int j = 0; j < node_cnt; ++j){
			node_t* node = nets.nets[i].vrtx.nodes[j];
			point_t shift = point_t_dif((*node).next, (*node).coord);
			st->mid_shift[cnt] = point_t_sum(st->mid_shift[cnt], shift);
			++cnt;
		}
	}
	++(st->it_cnt);
}

double statistic_t_get_max_mid_diviation(statistic_t st){
    double max_shift = 0;
    for (unsigned int l = 0; l < st.cnt_nodes; ++l){
        double len = SQR_LEN(st.mid_shift[l]);
        if (max_shift < len) max_shift = len;
    }
    return sqrt(max_shift)/st.it_cnt;
}

double statistic_t_get_full_mid_diviation(statistic_t st){
    double max_shift = 0;
    for (unsigned int l = 0; l < st.cnt_nodes; ++l)
        max_shift += SQR_LEN(st.mid_shift[l]);
    return sqrt(max_shift)/st.it_cnt;
}

void statistic_t_reset(statistic_t *st){
    memset(st->mid_shift, 0, st->cnt_nodes * sizeof(point_t));
    st->it_cnt = 0;
}

void set_initial_solving_params(world_t* world){
    collision_data_t_update(world->union_nets, &world->collision);
    --world->collision.counter;
    compute_nexts(world->union_nets, world->conditions.P, world->solver_data.delta, \
                        collision_t_get_state(world->collision), world->collision.box);
    update_statistic_t(world->dynamic_nets, &world->statistic);

    world->statistic.init_div = statistic_t_get_full_mid_diviation(world->statistic);
    statistic_t_reset(&world->statistic);
}

long double compute_dif_ms(struct timeval end, struct timeval start){
    long double seconds  = end.tv_sec  - start.tv_sec;
    long double useconds = end.tv_usec - start.tv_usec;
    long double mtime = (seconds) * 1000 + useconds / 1000.0;
    return mtime;
}

//eps пока никак не используется
void compute_nets_time(long double compute_time, world_t* world, int max_its)
{
    struct timeval start, end;
	gettimeofday(&start, NULL);
	int i = 0, crush = 0;
	for (i = 0; i < max_its; ++i){
        collision_data_t_update(world->union_nets, &world->collision);
        crush += compute_nexts(world->union_nets, world->conditions.P, world->solver_data.delta, \
                        collision_t_get_state(world->collision), world->collision.box);
        update_statistic_t(world->dynamic_nets, &world->statistic);
        update_nets(world->dynamic_nets);

        gettimeofday(&end, NULL);
		if (compute_dif_ms(end, start) >= compute_time - 0.8) break;
	}
	double res = statistic_t_get_full_mid_diviation(world->statistic);
	double rms = res / world->statistic.cnt_nodes;
//	static int n = 0;
//	n++;
//	if (!(n % 100))
//        printf("RMS shift per iter of node = %e mm, relation = %lg, cr = %d / %d\n", rms , res / world->statistic.init_div, crush, i+1);
	statistic_t_reset(&world->statistic);

}

point_t get_init_shift(point_t points[3]){			//не возвращает ошибки
	point_t normal = OR_AREA(points[0], points[1], points[2]);
	if (EQ_ERR(normal, ZERO(), 10 * DBL_EPSILON))
		{printf("ERROR: get_init_shift(...)\n"); return get_zero_point();}
	NORM_S(&normal);
	if (DOT(points[0], normal) <= 10 * DBL_EPSILON) return normal;
	return get_zero_point();
}

void deform_net(net_t net, point_t init_points[3], point_t final_points[3]){
	point_t shift = get_init_shift(init_points);
	point_t new_init_points[3];
	for (int i = 0; i < 3; i++)
		new_init_points[i] = point_t_sum(shift, init_points[i]);
	for (unsigned int i = 0; i < net.vrtx.count; i++)
		net.vrtx.nodes[i][0].next = point_t_sum(net.vrtx.nodes[i][0].coord, shift);

	matrix3D_t transit_from = matrix3D_t_get_from_points(new_init_points);
	matrix3D_t transit_to = matrix3D_t_get_from_points(final_points);
	matrix3D_t converse = matrix3D_t_mul(transit_to, matrix3D_t_inverse(transit_from));
	for (unsigned int i = 0; i < net.vrtx.count; i++)
		net.vrtx.nodes[i][0].next = matrix3D_t_vec_mul(converse, net.vrtx.nodes[i][0].next);
	update_net(net);
}

point_t find_line_to_net_intersection(net_t net, point_t line[2], double max_sqr_cntr_mass_dist, int* id){
	int e_cnt = net.elems.count;
	point_t inter[2] = {0};
	int flag = 0;
	//double res = 500000; int j = 0;
	for(int i = 0; i < e_cnt; i++){
		point_t cntr = net.elems.elems[i][0].cntr_mass;
		point_t projection = point_to_line_projection(cntr, line);
		//{point_t_dump(projection); point_t_dump(cntr); printf("dist = %lg\n", point_t_sqr_len(point_t_dif(cntr, projection)));}
		//{int cur = point_t_sqr_len(point_t_dif(cntr, projection)); if (cur < res) {res = cur; j = i;}}
		if (SQR_LEN(DIF(cntr, projection)) < max_sqr_cntr_mass_dist){
			//printf("HEY\n");
			point_t triangle[3];
			for (int j = 0; j < 3; j++)
				triangle[j] = net.elems.elems[i][0].vrts[j][0].coord;
			if (line_to_triangle_intersection(line, triangle, inter + (flag > 0))) {flag++;
				if (id) *id = i;
			}
			if (flag > 1 && DOT(DIF(line[1], line[0]), DIF(inter[1], inter[0])) > 0) {
				inter[0] = inter[1];
			}
		}
	}
	//printf("flag = %d\n", flag);
	//printf("res = %lg, max = %lg, i = %d\n", res, max_sqr_cntr_mass_dist, j);
	//if (flag == 2 && point_t_scal_mul(point_t_dif(line[1], line[0]), point_t_dif(inter[1], inter[0])) > 0) return inter[1];
	return inter[0];
}

double get_max_sqr_cntr_mass_dist(net_t net)
{
	double d = 0;
	unsigned int n_elems = net.elems.count, j;
	for (j = 0; j < n_elems; j++){
		elem_t* elem = net.elems.elems[j];
		point_t cntr = (*elem).cntr_mass;
		double x1 = SQR_LEN(DIF(cntr, (*((*elem).vrts[0])).coord));
		double x2 = SQR_LEN(DIF(cntr, (*((*elem).vrts[1])).coord));
		double x3 = SQR_LEN(DIF(cntr, (*((*elem).vrts[2])).coord));
		double cur_d = max_3d(x1, x2, x3);
		if (cur_d > d) d = cur_d;
	}
	return d;
}

void flat_obj_projection_to_net(net_t flat, net_t net){		//нужно поменять алгоритм выбора нормали
	int n_cnt = flat.vrtx.count;
	point_t flat_normal = flat.elems.elems[0][0].or_area;
	SCAL_S(-1.0 / LEN(flat_normal), &flat_normal);
	double max_dist = get_max_sqr_cntr_mass_dist(net);
	for (int i = 0; i < n_cnt; i++)
		if (is_fix(flat.vrtx.nodes[i][0].state)){
			node_t* node = flat.vrtx.nodes[i];
			point_t line[2];
			line[0] = node[0].coord;
			line[1] = SUM(line[0], flat_normal);
			node[0].next = find_line_to_net_intersection(net, line, max_dist, NULL);
		}
	update_net(flat);
}

//the following function works correctly only for free edges,
// that have only alone way, connecting all nodes of the edge
double net_t_get_len_of_way(net_t net, int (*is_node_belongs_way)(unsigned int)){
	node_t** nodes = net.vrtx.nodes;
	unsigned int n_cnt = net.vrtx.count;
	data_t way = data_t_construct(1, (int)(n_cnt * 0.05) + 1);
	for(unsigned int i = 0; i < n_cnt; i++)
		if ((*is_node_belongs_way)((*(nodes[i])).state)) add_elem_to_data(way, i, 0);

	double len = 0;
	int n_way = way.busy_len[0];
	for (int i = 0; i < n_way; i++){
		unsigned int n_spr = (*(nodes[way.data[0][i]])).cnt_springs;
		for (unsigned int j = 0; j < n_spr; j++){
			int m = 0;
			for (int k = 0; (k < n_way) && (m < 2); k++){
				if (k == i) continue;
				unsigned int sp_id = (*(nodes[way.data[0][i]])).springs_id[j];
				spring_t* spring = net.springs.springs[sp_id];
				if (spring_t_node_belong(spring, nodes[way.data[0][k]])){
					double spr_len = (*spring).l;
					len += spr_len;
					m++;
				}
			}
		}
	}
	data_t_destruct(way);

	return len / 2;
}

double net_t_get_len_free_edge(net_t net){
	return net_t_get_len_of_way(net, is_free_edge);
}

double net_t_get_len_fix_edge(net_t net){
	return net_t_get_len_of_way(net, is_fix);
}

double new_fix_len(net_t leaflet, point_t init_points[3], net_t aorta, point_t final_points[3], int flag){
	net_t newleaf = leaflet;
	if (!flag) newleaf= cpy_net(leaflet);
	deform_net(newleaf, init_points, final_points);
	flat_obj_projection_to_net(newleaf, aorta);
	double res = net_t_get_len_fix_edge(newleaf);
	if (!flag) net_t_destruct(newleaf);
	return res;
}

point_t third_point(point_t p1, point_t p2, point_t shift){
	return ADD(shift, 0.5, SUM(p1, p2));
}

line_t arrange_pnts_on_line(line_t line, int count){
	line_t res = {};
	if (count < 2) return perror("Wrong count"), res;
	double len = 0;
	for (int i = 0; i < line.pnt_cnt - 1; i++)
		len += LEN(DIF(line.line[i + 1], line.line[i]));

	res.pnt_cnt = count;
	res.line = (point_t*)calloc(count, sizeof(point_t));
	int cur = 0;
	double cur_len[2] = {0, 0};
	res.line[0] = line.line[0];
	for (int i = 1; i < count; i++){
		double place = len / (count - 1) * i;
		while (cur_len[1] < place){
			cur_len[0] = cur_len[1];
			cur_len[1] += LEN(DIF(line.line[cur + 1], line.line[cur]));
			cur++;
		}
		double coef = (place - cur_len[0]) / (cur_len[1] - cur_len[0]);
		res.line[i] = SCAL(coef, DIF(line.line[cur], line.line[cur - 1]));
		res.line[i] = SUM(line.line[cur - 1], res.line[i]);
	}

	return res;
}

int get_next_fix_node(net_t net, unsigned int cur_id, unsigned int prev_id){
	node_t* node = net.vrtx.nodes[cur_id];
	for (unsigned int i = 0; i < node[0].cnt_springs; i++){
		spring_t* spr = net.springs.springs[node[0].springs_id[i]];
		node_t* neigh = spring_t_get_other_end(spr, node);
		if (is_fix(neigh[0].state) && neigh[0].id != prev_id)
			return neigh[0].id;
	}
	return -1;
}

data_t get_fix_bnd(net_t net, unsigned int start_node_id){
	data_t fix_bnd = data_t_construct(1, 100);
	int cur_id = start_node_id;
	int prev_id = cur_id;
	while (cur_id > 0){
		add_elem_to_data(fix_bnd, (unsigned int)cur_id, 0);
		int new_id = get_next_fix_node(net, (unsigned int)cur_id, (unsigned int)prev_id);
		prev_id = cur_id;
		cur_id = new_id;
	}
	return fix_bnd;
}

void sew_leaflet_to_bnd(net_t leaflet, unsigned int start_node_id, line_t bnd){
	data_t fix_bnd = get_fix_bnd(leaflet, start_node_id);
	line_t newcoord = arrange_pnts_on_line(bnd, fix_bnd.busy_len[0]);
	for (int i = 0; i < fix_bnd.busy_len[0]; i++)
		leaflet.vrtx.nodes[fix_bnd.data[0][i]][0].next = newcoord.line[i];
	update_net(leaflet);
	data_t_destruct(fix_bnd);
	free(newcoord.line);
}

/*void sew_leaflet_to_aorta(	net_t leaflet, point_t init_attachment[2], point_t init_direction,\
							net_t aorta,   point_t attachment[2], 	   point_t blood_direction){
	double len = net_t_get_len_fix_edge(leaflet);
	printf("init_fix_len = %lg\n", len);
	point_t init_points[3], final_points[3];
	init_points[0] = init_attachment[0], init_points[1] = init_attachment[1];
	final_points[0] = attachment[0], final_points[1] = attachment[1];
	init_points[2] = third_point(init_points[0], init_points[1], init_direction);
	double coef = 1;
	double coef_p = -1;
	double coef_m = -1;
	double err = 1;
	const double constraint = 0.01;
	int flag = 1;
	while (flag){
		flag = (fabs(err) > constraint);
		point_t newdirect = blood_direction;
		point_t_coef_mul(coef, &newdirect);
		final_points[2] = third_point(final_points[0], final_points[1], newdirect);
		double res = new_fix_len(leaflet, init_points, aorta, final_points, (fabs(err) < constraint));
		printf("cur_fix_len = %lg\n", res);
		err = (res - len) / len;
		double oldcoef = coef;
		if (err > 0 && (coef_p > coef || coef_p < 0)) coef_p = coef;
		if (err < 0 && coef_m < coef) coef_m = coef;
		if (coef_p > 0 && coef_m > 0) coef = (coef_p + coef_m) / 2;
		else coef = coef  / (err + 1);
		if ((fabs(err) < constraint)) coef = oldcoef;
	}

}*/

void sew_leaflet_to_aorta(	net_t leaflet, point_t init_attachment[2], point_t init_direction,\
							net_t aorta,   point_t attachment[2], 	   point_t blood_direction,\
							unsigned int start_node_id, line_t bnd){
	point_t init_points[3], final_points[3];
	init_points[0] = init_attachment[0], init_points[1] = init_attachment[1];
	final_points[0] = attachment[0], final_points[1] = attachment[1];
	init_points[2] = third_point(init_points[0], init_points[1], init_direction);
	final_points[2] = blood_direction;
	normalize(&final_points[2]);
	point_t_coef_mul(point_t_length(init_direction), &final_points[2]);
	final_points[2] = third_point(final_points[0], final_points[1], final_points[2]);
	deform_net(leaflet, init_points, final_points);

	sew_leaflet_to_bnd(leaflet, start_node_id, bnd);
}
//############computation###############################################

