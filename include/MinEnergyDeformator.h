#ifndef MINENERGYDEFORMATOR_H
#define MINENERGYDEFORMATOR_H

#include <cstring>
#include <iostream>
#include "nets.h"
#include "precomputation.h"
#include "computation.h"
#include "save-data.h"
#include <stdlib.h>
#include <math.h>
#include <stdio.h>
#include <sstream>
#include <vector>
#include <set>
#include <algorithm>
#include <gsl/gsl_multimin.h>
#include <chrono>
#include <memory>

//нужно умными указателями сделать сохранение функтора

class MinEnergyDeformator{
private:
    net_t& m_mesh;
public:
    class EnergyFunctor{
        static inline double zero(const net_t& mesh, int node_id, double* df) { return 0; }
    private:
        using Energy = double(*)(const net_t& mesh, int node_id, double* df);
        Energy func;
    public:
        EnergyFunctor( double(*func)(const net_t& mesh, int node_id, double* df) = zero): func{func} {};
        virtual double operator()(const net_t& mesh, int node_id, double* df)
        { return func(mesh, node_id, df); }
    };

    // this functions should save gradient into df, if df != nullptr
    using NodeEnergy =  std::shared_ptr<EnergyFunctor>; //df[3]
    using EdgeEnergy = std::shared_ptr<EnergyFunctor>; //df[6]
    using DigedralEnergy = std::shared_ptr<EnergyFunctor>; //df[12]

private:
    std::vector<NodeEnergy> m_node_energy;
    std::vector<EdgeEnergy> m_edge_energy;
    std::vector<DigedralEnergy> m_digedral_energy;

    std::vector<int> m_node_map;
    int m_N_vars = 0;
    //std::vector<std::array<point_t, 2>> m_digedrals;

public:
    MinEnergyDeformator(net_t& mesh): m_mesh{mesh}
    {
        set_node_map();
    }

    void addNodeEnergyComponent(NodeEnergy f) { m_node_energy.push_back(f); }
    void addEdgeEnergyComponent(EdgeEnergy f) { m_edge_energy.push_back(f); }
    void addDigedralEnergyComponent(DigedralEnergy f) { m_digedral_energy.push_back(f); }
    double getEnergy() { return compute(); }
private:

    inline void set_coords_to_mesh(const gsl_vector* X){
        for (int j = 0, i = 0; j < m_mesh.vrtx.count; ++j){
            if (is_fix(m_mesh.vrtx.nodes[j]->state)){
                continue;
            }
            for (int k = 0; k < 3; ++k)
                m_mesh.vrtx.nodes[j]->coord.coord[k] = gsl_vector_get(X, i++);
        }
    }

    inline void set_node_map(){
        int i = 0;
        m_node_map.resize(m_mesh.vrtx.count);
        for (int j = 0; j < m_mesh.vrtx.count; ++j){
            if (is_fix(m_mesh.vrtx.nodes[j]->state)){
                m_node_map[j] = -1;
                continue;
            }
            m_node_map[j] = i;
            i++;
        }
        m_N_vars = 3 * i;
    }

    double compute(){
        double f = 0;
        for (int i = 0; i < m_mesh.springs.count; ++i){
            spring_t& sp = *m_mesh.springs.springs[i];
            if (is_fix(sp.ends[0]->state) && is_fix(sp.ends[1]->state)) continue;
            for (auto& en: m_edge_energy)
                f += en->operator()(m_mesh, i, nullptr);
            if (sp.isdigedral)
                for (auto& en: m_digedral_energy)
                    f += en->operator()(m_mesh, i, nullptr);
        }

        for (int i = 0; i < m_mesh.vrtx.count; ++i)
            for (auto& en: m_node_energy)
                f += en->operator()(m_mesh, i, nullptr);
        return f;
    }

    inline void set_templated(int* id, int n, gsl_vector *df, double* ddf){
        for (int i = 0; i < n; ++i)
            if (id[i] >= 0)
                for (int j = 0; j < 3; ++j)
                    gsl_vector_set(df, id[i] * 3 + j,gsl_vector_get(df, id[i] * 3 + j) + ddf[i*3 + j]);
    }

    void compute_df(gsl_vector *df){
        gsl_vector_set_all(df, 0);

        for (int i = 0; i < m_mesh.springs.count; ++i) {
            spring_t& sp = *m_mesh.springs.springs[i];
            {
                double loc_df[6] = {};
                for (auto &en: m_edge_energy)
                    en->operator()(m_mesh, i, loc_df);
                int id[2] = {m_node_map[sp.ends[0]->id], m_node_map[sp.ends[1]->id]};
                set_templated(id, 2, df, loc_df);
            }

            if (sp.isdigedral) {
                double loc_df[12] = {};
                for (auto& en: m_digedral_energy)
                    en->operator()(m_mesh, i, loc_df);
                int id[4] = {m_node_map[sp.dihedral[0]->id], m_node_map[sp.dihedral[1]->id],
                             m_node_map[sp.ends[0]->id],     m_node_map[sp.ends[1]->id]};
                set_templated(id, 4, df, loc_df);
            }
        }

        for (int i = 0; i < m_mesh.vrtx.count; ++i){
            double loc_df[3] = { 0 };
            for (auto& en: m_node_energy)
                en->operator()(m_mesh, i, loc_df);
            set_templated(&m_node_map[m_mesh.vrtx.nodes[i]->id], 1, df, loc_df);
        }
    }

    static void _fdf4gsl(const gsl_vector *x, void *params, double *f, gsl_vector *df)
    {
        MinEnergyDeformator* m = (MinEnergyDeformator*)params;
        m->set_coords_to_mesh(x);
        *f = m->compute();
        m->compute_df(df);
    }

    static double _f4gsl (const gsl_vector *x, void *params)
    {
        MinEnergyDeformator* m = (MinEnergyDeformator*)params;
        m->set_coords_to_mesh(x);
        return m->compute();
    }

    static void _df4gsl (const gsl_vector *x, void *params, gsl_vector *df)
    {
        MinEnergyDeformator* m = (MinEnergyDeformator*)params;
        m->set_coords_to_mesh(x);
        m->compute_df(df);
    }

public:
    int find_minimum_energy_df(int freq = 1, double step_sz = 1.e-4, double tol = 1.e-4, double epsabs = 1.e-2, int maxits = 50000, double time/*secs*/ = 150){
        size_t iter = 0;
        int status;

        const gsl_multimin_fdfminimizer_type *T;
        gsl_multimin_fdfminimizer *s;

        void* par = this;

        int g_N = m_N_vars;
        gsl_vector *x = gsl_vector_alloc (g_N);
        for (int j = 0, i = 0; j < m_mesh.vrtx.count; ++j){
            if (is_fix(m_mesh.vrtx.nodes[j]->state)){
                continue;
            }
            for (int k = 0; k < 3; ++k)
                gsl_vector_set(x, i++, m_mesh.vrtx.nodes[j]->coord.coord[k]);// + (k == 2) ? 1.0 : 0.0);
        }

        gsl_multimin_function_fdf my_func;

        my_func.n = g_N;
        my_func.f = _f4gsl;
        my_func.df = _df4gsl;
        my_func.fdf = _fdf4gsl;
        my_func.params = par;

        T = gsl_multimin_fdfminimizer_conjugate_fr;
        s = gsl_multimin_fdfminimizer_alloc (T, g_N);

        gsl_multimin_fdfminimizer_set (s, &my_func, x, step_sz, tol);

        class _Timer
        {
        private:
            using clock_t = std::chrono::high_resolution_clock;
            using second_t = std::chrono::duration<double, std::ratio<1> >;
            std::chrono::time_point<clock_t> m_beg;
        public:
            _Timer() : m_beg(clock_t::now()){}
            void reset() { m_beg = clock_t::now(); }

            double elapsed() const
            {
                return std::chrono::duration_cast<second_t>(clock_t::now() - m_beg).count();
            }
        };

        _Timer t;
        do
        {
            iter++;
            status = gsl_multimin_fdfminimizer_iterate (s);

            if (status)
                break;

            status = gsl_multimin_test_gradient (s->gradient, epsabs);

            if (status == GSL_SUCCESS)
                printf ("Minimum found at: %5lu: f() = %10.5f\n", iter, s->f);

            if ((freq > 0) && !(iter % freq))
                printf ("%5lu: f() = %10.5f\n", iter, s->f);

        }
        while (status == GSL_CONTINUE && iter < maxits && t.elapsed() < time);
        if (status != GSL_SUCCESS)
            printf ("Success didn't reached status = %d: %5lu: f() = %10.5f\n", status, iter, s->f);

        gsl_multimin_fdfminimizer_free (s);
        gsl_vector_free (x);

        return status;
    }

};

void set_default_length_constr(MinEnergyDeformator& m, double weight);
void set_default_digedral_angle_constr(MinEnergyDeformator& m, double weight, double scale = 1.0, bool convexity = true);
void set_plane_constr(MinEnergyDeformator& m, double weight, const point_t& normal, double b);
void set_isotrop_force(MinEnergyDeformator& m, double weight, const point_t& force);

#endif