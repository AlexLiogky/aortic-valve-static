#include "elastic_forces.h"

static point_t elastic_force_REINFORCING_(net_t net, node_t* node, double* prms){
    double mu = prms[0] / 3.0;
    double Jm = prms[1];

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

        point_t ff = ZERO();
        //assert(LEN(ff) < 10 * DBL_EPSILON);
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

int set_elastic_force_REINFORCING(net_t net, double delta, void* prms){
    //return set_elastic_force_dummy(net, delta, prms, elastic_force_REINFORCING_);
    double* to = get_params(prms);
    double mu = to[0] / 3;
    double Jm = to[1];
    for (int i = 0; i < net.elems.count; ++i){
        elem_t* e = net.elems.elems[i];
        double Ap = e->area_0;
        double Aq = LEN(e->or_area);
        double J = Aq / Ap;
        node_t** n = e->vrts;
        double* m = to + 2 + 9 * i;
        double trC = 0;
        for (int j = 0; j < 3; ++j)
            for (int k = 0; k < j; ++k)
                trC += 2 * DOT(n[j]->coord, n[k]->coord) * m[j*3 + k];
        for (int j = 0; j < 3; ++j)
            trC += SQR_LEN(n[j]->coord) * m[j * 3 + j];
        double _scale = mu / (1 - (trC + 1 / J / J - 3) / Jm) * delta;
        for (int j = 0; j < 3; ++j){
            double scale = _scale * n[j]->h;
            for (int jj = 0; jj < 3; ++jj)
                ADD_S(&(n[j]->next), -scale * Ap * m[j * 3 + jj], n[jj]->coord);
            ADD_S(&(n[j]->next), scale / 2 / J / J / J, get_ortho_vector(n[(j + 1)%3]->coord, n[(j + 2)%3]->coord, n[j]->coord));
        }
    }
}


int precomp_elastic_force_REINFORCING(net_t net, double* prms, void** data){
    //return dummy_saving_prms(data, prms, 2);
    int cnt = 2 + 9 * net.elems.count;
    double* to = resize(data, cnt);
    int off = 0;
    to[off++] = prms[0];
    to[off++] = prms[1];
    double mu = prms[0] / 3;
    for (int i = 0; i < net.elems.count; ++i){
        elem_t* e = net.elems.elems[i];
        node_t** n = e->vrts;
        double Ap = e->area_0;
        point_t D[3] = {};
        for (int j = 0; j < 3; ++j)
            D[j] = SCAL(1.0 / 2.0 / Ap, get_ortho_vector(n[(j + 1)%3]->initial, n[(j + 2)%3]->initial, n[(j)%3]->initial));
        double m[3][3];
        for (int ii = 0; ii < 3; ++ii)
            for (int jj = ii + 1; jj < 3; ++jj){
                m[ii][jj] = DOT(D[ii], D[jj]);
                m[jj][ii] = m[ii][jj];
            }
        for (int ii = 0; ii < 3; ++ii)
            m[ii][ii] = SQR_LEN(D[ii]);
        for (int ii = 0; ii < 9; ++ii)
            to[off++] = m[ii/3][ii%3];
    }

    return 0;
}