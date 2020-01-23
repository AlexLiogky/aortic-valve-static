#include "elastic_forces.h"

static point_t elastic_force_NEOGOOK_(net_t net, node_t* node, double* prms){
    double mu = prms[0] / 3;

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
        //assert(flag);

        double Ap = triag->area_0;
        point_t D[3] = {};
        for (int j = 0; j < 3; ++j)
            D[j] = SCAL(1.0 / 2.0 / Ap, get_ortho_vector(n[(j + 1)%3]->initial, n[(j + 2)%3]->initial, n[(j)%3]->initial));
        double Aq = LEN(triag->or_area);

        double J = Aq / Ap;
        for (int j = 0; j < 3; ++j)
            ADD_S(&f, mu * Ap * DOT(D[0], D[j]), n[j]->coord);

        ADD_S(&f, -mu / 2 / J / J / J, get_ortho_vector(n[1]->coord, n[2]->coord, n[0]->coord));
    }

    SCAL_S(-node->h, &f);

    return f;

}

int set_elastic_force_NEOGOOK(net_t net, double delta, void* prms){
    return set_elastic_force_dummy(net, delta, prms, elastic_force_NEOGOOK_);
}


int precomp_elastic_force_NEOGOOK(net_t net, double* prms, void** data){
    return dummy_saving_prms(data, prms, 1);
}