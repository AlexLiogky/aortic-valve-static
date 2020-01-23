#include "elastic_forces.h"
#include <math.h>

static point_t elastic_force_TRQS_(net_t net, node_t* node, double* prms){     //St Venan
    double E = prms[0];      //Young modulus, Pa
    double nu = prms[1];           //Poisson ratio

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
            //assert(spr != NULL);
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

    return f;   //10^-6 * H
}

int set_elastic_force_TRQS(net_t net, double delta, void* prms){
    return set_elastic_force_dummy(net, delta, prms, elastic_force_TRQS_);
}


int precomp_elastic_force_TRQS(net_t net, double* prms, void** data){
    return dummy_saving_prms(data, prms, 2);
}