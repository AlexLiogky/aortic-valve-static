#include "elastic_forces.h"
#include "precomputation.h"

static double spring_t_set_coef_(net_t net, spring_t* obj){
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

int precomp_elastic_force_MSM(net_t net, double* prms, void** data){
    const double _E_0 = prms[0];
    const double _E_1 = prms[1];
    const double _Eps_lim = prms[2];
    double* to = resize(data, 3 * net.springs.count);
    to[0] = _Eps_lim;
    for (int i = 0, cnt = net.springs.count; i < cnt; ++i){
        spring_t* spr = net.springs.springs[i];
        double coef = spring_t_set_coef_(net, spr);
        double save[2] = {_E_0 * coef, _E_1 * coef};
        to[2*i + 1] = save[0];
        to[2*i + 2] = save[1];
    }

    return 0;
}

static double stress_t_get_force_(double deformation, double* precomp){
    double force = 0;
    double eps_lim = precomp[0];
    double sigma1 = precomp[1];
    double sigma2 = precomp[2];
    int sign = (deformation > 0) ? 1 : -1;
    force += (fabs(deformation) < eps_lim) ? sigma1 * deformation : sigma1 * eps_lim * sign;
    force += (fabs(deformation) > eps_lim) ? sigma2 * (deformation - sign * eps_lim) : 0;
    return force;
}

int set_elastic_force_MSM(net_t net, double delta, void* prms){
    double* to = get_params(prms);
    const double _Eps_lim = to[0];
    for (int i = 0, cnt = net.springs.count; i < cnt; ++i) {
        spring_t* spr = net.springs.springs[i];
        double l = spr->l, l_0 = spr->l_0;
        double eps = (l - l_0) / l_0;
        double param[3] = {_Eps_lim, to[2*i + 1], to[2*i + 2]};
        double force = stress_t_get_force_(eps, param);
        ADD_S(&(spr->ends[0]->next), force * delta, spr->direction);
        ADD_S(&(spr->ends[1]->next), -force * delta, spr->direction);
    }
}