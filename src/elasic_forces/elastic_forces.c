#include "elastic_forces.h"

double* resize(void** data, long nmembs){
    long long cur_sz = ((long long*)(*data))[0];
    long long predicted = 2 * sizeof(long long) + sizeof(double) * nmembs;
    if (predicted > cur_sz) {
        *data = realloc(*data, predicted);
        ((long long*)(*data))[0] = predicted;
    }
    double* to = get_params(*data);
    return to;
}

double* get_params(void* data){
    return (double*)(data + 2 * sizeof(long long));
}

int dummy_saving_prms(void** data, double* prms, int cnt){
    double* to = resize(data, cnt);
    for (int i = 0; i < cnt; ++i)
        to[i] = prms[i];
    return 0;
}

int set_elastic_force_dummy(net_t net, double delta, void* prms, point_t (*f)(net_t net, node_t* node, double* prms)){
    double* to = get_params(prms);
    for (int i = 0, cnt = net.vrtx.count; i < cnt; ++i){
        node_t* n = net.vrtx.nodes[i];
        point_t ff = f(net, n, to);
        ADD_S(&n->next, delta, ff);
    }
    return 0;
}