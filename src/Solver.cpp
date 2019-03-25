#include "Solver.h"

Solver::Solver(nets_t& dynamic_nets, nets_t& static_nets, wrld_cnd_t& cond, solver_t& solver_data, double drop_thr, double max_shft):
    m_dynamic_nets{dynamic_nets},
    m_static_nets{static_nets},
    m_union_nets{create_union_net(dynamic_nets, static_nets)},
    m_conditions{cond},
    m_solver_data{solver_data},
    m_statistic{statistical_data_construct(dynamic_nets)},
    m_allow_shift{drop_thr},
    m_max_shift{max_shft}
{
    m_collision.reserve(m_union_nets.count);
    for (int i = 0; i < m_union_nets.count; ++i)
        m_collision.push_back(new Net_Wraper(m_union_nets.nets[i], m_conditions.P, m_solver_data.delta));

    set_initial_solving_params();
}

void Solver::set_initial_solving_params()
{
    predictMotion();
    findCollisions();
    solveCollisions();
    update_statistic_t(m_dynamic_nets, &m_statistic);

    m_statistic.init_div = statistic_t_get_full_mid_diviation(m_statistic);
    statistic_t_reset(&m_statistic);
}

double Solver::predictMotion(){
    double shift = 0;
    for (auto& net: m_collision){
        double max_net_shift = net->computeFreeNexts();
        if (shift < max_net_shift) shift = max_net_shift;
    }
    char ret = 0;
    if (shift > m_max_shift) {
		relaxation(m_union_nets, sqrt(m_allow_shift / shift));
		ret++;
	}
    for (auto& net: m_collision){
        net->updateCollisionInfo();
    }
    return ret;
}

void Solver::findCollisions(){
    for (int i = 1; i < m_collision.size(); ++i)
        for (int j = 0; j < i; ++j)
        {
            m_collision[i]->CollisionHandler(m_collision[j]);
        }
}

void Solver::solveCollisions(){
    for (auto& body: m_collision)
        Net_Wraper::solveSContacts(body);
}

void Solver::compute_nets_time(long double compute_time, int max_its)
{
    struct timeval start, end;
	gettimeofday(&start, NULL);
	int i = 0, crush = 0;
	for (i = 0; i < max_its; ++i){
        //collision_data_t_update(m_union_nets, &world->collision);
        crush += predictMotion();
        findCollisions();
        solveCollisions();

        update_statistic_t(m_dynamic_nets, &m_statistic);
        update_nets(m_dynamic_nets);

        gettimeofday(&end, NULL);
		if (compute_dif_ms(end, start) >= compute_time - 0.8) break;
	}
	double res = statistic_t_get_full_mid_diviation(m_statistic);
	double rms = res / m_statistic.cnt_nodes;
	//printf("RMS shift per iter of node = %e mm, relation = %lg, cr = %d / %d\n", rms , res / world->statistic.init_div, crush, i+1);
	statistic_t_reset(&m_statistic);
}


