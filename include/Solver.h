#ifndef SOLVER_H
#define SOLVER_H
#include "computation.h"
#include "nets.h"
#include <vector>
#include "Net_Wraper.h"

using namespace std;

class Solver
{
public:
        nets_t& getDynamicNets() { return m_dynamic_nets; }
        nets_t& getStaticNets() { return m_static_nets; }
        Solver(nets_t& dynamic_nets, nets_t& static_nets, wrld_cnd_t& cond, solver_t& solver_data, double drop_thr, double max_shft);
        void compute_nets_time(long double compute_time = 1000.0/60, int max_its = 10000);

private:
    void set_initial_solving_params();
protected:
    double predictMotion();
    void findCollisions();
    void solveCollisions();

private:
    nets_t m_dynamic_nets;
    nets_t m_static_nets;
    nets_t m_union_nets;
    std::vector<Net_Wraper*> m_collision;
    wrld_cnd_t m_conditions;
    solver_t m_solver_data;
    statistic_t m_statistic;
	double m_allow_shift = 1;
	double m_max_shift;

};

#endif // SOLVER_H
