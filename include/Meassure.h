#ifndef MEASSURE_H_INCLUDED
#define MEASSURE_H_INCLUDED

#include "nets.h"
#include "save-data.h"
#include "format_in.h"
#include "precomputation.h"
#include "World.h"
#include "InputProcessor.h"
#include <string>
#include <sstream>
#include <fstream>
#include <iostream>
#include <sstream>
#include <map>
#include <set>
#include <algorithm>
#include <iterator>
#include <cmath>

#include <gsl/gsl_multimin.h>
#include <chrono>

#include <CGAL/Simple_cartesian.h>
#include <CGAL/linear_least_squares_fitting_3.h>


using namespace std;

//  возможно в будущем подобные функции будут использоваться в real time debug версии, поэтому потом имеет смысл сделать
// эффективную по памяти и алгоритму версию этого класса,
// в первую очередь можно модифицировать getCollision чтобы он возвращал более удобные данные
class CoaptHeight{
private:
    struct ComparePoints{
        bool operator()(const point_t& p1, const point_t& p2) const
        {
            for (int i = 0; i < 3; ++i)
                if (p1.coord[i] > p2.coord[i]) return false;
                else if (p1.coord[i] < p2.coord[i]) return true;
            return false;
        }
    };
    struct StrongCompareNodes{
        bool operator()(const node_t* p1, const node_t* p2) const
        {
            return ComparePoints()(p1->coord, p2->coord);
        }
    };
    struct CompareNodes{
        bool operator()(const node_t* p1, const node_t* p2) const
        {
            return p1->id < p2->id;
        }
        bool operator()(const node_t& p1, const node_t& p2) const
        {
            return p1.id < p2.id;
        }
    };

    struct RefineElem{
        auto get_other_node(elem_t* elem, node_t* n1, node_t* n2);
        auto get_neigh_elem(const net_t net, const elem_t* elem, const node_t* n1, const node_t* n2);
        pair<array<elem_t*, 3>, int> compute_shared_elems(net_t net, elem_t* e);
        int get_state(int state1, int state2);
        elem_t* e_construct(node_t* node1, node_t* node2, node_t* node3, unsigned int id, elem_t* on);
        spring_t* spr_construct(node_t* node1, node_t* node2, unsigned int id);
        template<int N>
        void set_data_to_existing_node(node_t* n, set<int>& excluding, set<int> what) {
            if constexpr (N == 0)
                set_data_to_existing_node_0(n, excluding, what);
            else if constexpr (N == 1)
                set_data_to_existing_node_1(n, excluding, what);
        }
        void set_data_to_existing_node_0(node_t* n, set<int>& excluding, set<int> what);
        void set_data_to_existing_node_1(node_t* n, set<int>& excluding, set<int> what);
        void compute_new_indeces();
        void create_new_objects();
        void update_old_nodes();
        void insert_new_objs_into_net();
        void prepare_elems();
        void clear();
        void test();
        vector<elem_t*> refine_elem(net_t& net, elem_t* e);

        net_t* m_net;
        elem_t* m_e;
        bool filled = false;
        map<int, pair<elem_t*, node_t*>> m_neigh_elems;
        vector<int> m_nodes_id;
        vector<int> m_springs_id;
        vector<int> m_elems_id;
        vector<node_t*> m_new_nodes;
        vector<spring_t*> m_new_sprs;
        vector<elem_t*> m_new_elems;
        map<int, set<int>> m_nds_els_id;
        map<int, set<int>> m_nds_sprs_id;
    };

    struct To_colData{
        solver_t solver_data;
        wrld_cnd_t conditions;
        double colission_margin;
        double computive_mrg;
    };

    typedef double                      FT;
    typedef CGAL::Simple_cartesian<FT>  K;
    typedef K::Plane_3                  Plane;
    typedef K::Point_3                  Point_3;
    typedef K::Vector_3                 Vector_3;
    typedef K::Vector_2                 Vector_2;
    typedef K::Point_2                  Point_2;
    typedef K::Segment_2                Segment_2;
    typedef K::Line_2                   Line_2;
    typedef K::Triangle_3               Triangle;

private:
    void set_colData();
    template<typename ElemIterator, typename ElemToPoints>
    void _savetoSTL(std::ofstream& stl, ElemIterator& eit, const ElemToPoints& e_to_pnts, string name = ""){
        stl << "solid " << name << std::endl;
        for (auto& j: eit){
            auto pnts = e_to_pnts(j);
            point_t n = NORM(OR_AREA(pnts[0], pnts[1],pnts[2]));
            stl << "facet normal " << n.coord[0] << " " << n.coord[1] << " " << n.coord[2] << endl;
            stl << "outer loop" << endl;
            for (int k = 0; k < 3; ++k) {
                double* p = pnts[k].coord;
                stl << "vertex " << p[0] << " " << p[1] << " " << p[2] << endl;
            }
            stl << "endloop" << endl;
            stl << "endfacet" << endl;
        }
        stl << "endsolid " << name << std::endl;
    }
    template<typename NodeContainer, typename NodeToPoints>
    void _savetoOFF(std::ofstream& off, const NodeContainer& nct, const NodeToPoints& n_to_dbl_vec){
        off << "OFF\n\n";
        off << nct.size() << " " << "0" << " 0\n";
        for (auto& j: nct)
        {
            auto pnt = n_to_dbl_vec(j);
            for (auto k: pnt)
                off << k << " ";
            for (int k = pnt.size(); k < 3; ++k)
                off << 0 << " ";
            off << "\n";
        }
    }
    template <int DD, typename T1, typename T2, typename Convert>
    auto get_linear_least_squares_fitting_plane(T1& data1, T2& data2, Convert conv){
        auto make_set = [](const T1& t0, const auto& ...xs){
            set<typename T1::value_type> s;
            auto my_insert = [&s](auto param) mutable {
                for (auto& i: param)
                    s.insert(i);
            };
            my_insert(t0);
            return (void)initializer_list<int>{((void)my_insert(xs), 0)...}, s;
        };
        auto s = make_set(data1, data2);
        vector<decltype(conv(data1[0]))> objs; objs.reserve(s.size());
        for (auto& i: s)
            objs.push_back(conv(i));
        Plane plane;
        linear_least_squares_fitting_3(objs.begin(), objs.end(), plane, CGAL::Dimension_tag<DD>());
        return plane;
    }
    template<typename NodeConverter, typename CheckEdgeExist>
    vector<array<int, 3>> convert_edges_to_elems(NodeConverter node_conv, CheckEdgeExist is_here, int mesh_n){
        vector<array<int, 3>> res;
        for (int i = 0, cnt = mesh.nets[mesh_n].elems.count; i < cnt; ++i){
            auto e = mesh.nets[mesh_n].elems.elems[i];
            array<int, 3> loop;
            bool flag = true;
            for (int j = 0; j < 3; ++j) {
                loop[j] = node_conv(e->vrts[j]);
                flag *= (loop[j] >= 0);
            }
            if (!flag) continue;
            array<int, 3> loop_cp = loop;
            sort(loop.begin(), loop.end());
            const array<pair<int, int>, 3> finder = {pair(0, 1), {0, 2}, {1, 2}};
            for (int j = 0; j < 3; ++j)
                flag *= (is_here({loop[finder[j].first], loop[finder[j].second]}));
            if (!flag) continue;
            res.push_back(loop_cp);
        }
        return res;
    }
    void computeFaces();
    Point_3 point_t_to_Point(const point_t& p0) const ;
    map<pair<int, int>, Plane> find_mid_planes(const int dim = 2);
    int set_mid_planes(const int dim = 2);
    Point_2 p2plane(const point_t& p, const Plane& plane);
    int project_coaption_to_planes(int state = 0);
    int project_main_direction(int state = 0);
    pair<Point_2, Point_2> get_minmax_Point_set(vector<Point_2>& data, Vector_2 dir);
    void set_borders();
    void set_triangle_edges();
    void set_all_edges();
    void compute_boundary_edges();
    static array<double, 2> _clever_field_directed_diam(const vector<Point_2>& pnts, const set<pair<int, int>>& edges, const Line_2& line, const Vector_2& dir);
    static double _field_directed_diam(const vector<Point_2>& pnts, const set<pair<int, int>>& edges, const Line_2& line, const Vector_2& dir);
    double field_directed_diam(const pair<int, int>& id, const Line_2& line, const Vector_2& dir);
    void compute_billow_plate();
    void find_refined_elems();
    vector<elem_t*> refine_elem(int mesh_n, elem_t* e);
    set<int> test_node_in_springs(node_t* n, net_t net);
    void compute_billowing_pnts();
    set<node_t*> common_nodes(int mesh_n);
    set<spring_t*> common_edges(int mesh_n);
    set<node_t*> common_nodes_with_connection(int mesh_n);
    static vector<Point_2> align_to_dir(const vector<Point_2>& pnts, const Vector_2& dir);
    static vector<Point_2> dealign_from_dir(const vector<Point_2>& aligned, const Vector_2& dir);

    struct SewTwoFields{
        SewTwoFields(vector<Point_2>& p1, set<pair<int, int>>& link1,
                     vector<Point_2>& p2, set<pair<int, int>>& link2, map<int, int>& common_pnts_id):
                     m_p1{p1}, m_p2{p2}, m_l1{link1}, m_l2{link2}, m_cn{common_pnts_id}
        {}
    private:
        void set_data(gsl_vector *x);
        void set_new_points(gsl_vector* x);
        double compute(const gsl_vector *x);
        void compute_df(const gsl_vector *x, gsl_vector *df);
        double compute_fdf(const gsl_vector *x, gsl_vector *df);
        static void _fdf4gsl(const gsl_vector *x, void *params, double *f, gsl_vector *df);
        static double _f4gsl (const gsl_vector *x, void *params);
        static void _df4gsl (const gsl_vector *x, void *params, gsl_vector *df);
    public:
        int find_minimum_energy_df(int freq = 1, double step_sz = 1.e-4, double tol = 1.e-4, double epsabs = 1.e-2, int maxits = 50000, double time = 150);

    private:
        vector<Point_2> &m_p1, &m_p2;
        set<pair<int, int>> &m_l1, &m_l2;
        map<int, int>& m_cn;
        const double m_wy = 10;
        map<int, int> m_map1, m_map2;
    };

    static tuple<vector<Point_2>, set<pair<int, int>>, map<int, int>>
            _get_connected_pnts(const vector<Point_2>& pnts, const set<pair<int, int>>& connect);
    void compute_plane_leaf_coaptation_profile(int mesh_n);
    void compute_leaf_profiles();
public:


private:
    nets_t mesh;
    bool mesh_changed = true;
    bool m_billow_plate_state = false;
    bool m_billowing_state = false;
    array<point_t, 3> m_billow_plate;
    array<point_t, 3> m_billow_pnt;
    array<double, 3> m_billowing;

    map<pair<int, int>, set<node_t*>> colData;
    map<pair<int, int>, vector<node_t*>> m_nodes;
    map<pair<int, int>, vector<Point_2>> m_projection;
    map<pair<int, int>, set<pair<int, int>>> m_edges;
    map<pair<int, int>, set<pair<int, int>>> m_boundary_edges;
    map<pair<int, int>, K::Vector_2> proj_dir;
    map<pair<int, int>, vector<array<int, 3>>> m_faces;
    map<pair<int, int>, pair<Point_2, Point_2>> m_borders;
    vector<vector<elem_t*>> m_no_coapt_elems;
    vector<map<elem_t*, bool>> m_refined_elems;
    RefineElem m_refiner;
    vector<map<int, pair<int, int>>> remap = {{{0, {0, 1}}, {1, {0, 2}}},
                                              {{0, {1, 0}}, {1, {1, 2}}},
                                              {{0, {2, 0}}, {1, {2, 1}}}};
    map<pair<int, int>, int> mapre = {{{0, 1}, 0}, {{0, 2}, 1},
                                      {{1, 0}, 0}, {{1, 2}, 1},
                                      {{2, 0}, 0}, {{2, 1}, 1}};
    map<pair<int, int>, Plane> m_planes;
    To_colData _to_colData;
    point_t axis;

    map<int, tuple<vector<Point_2>, set<pair<int, int>>, map<int, int>, set<pair<int, int>>>> m_leaflet_coapt_fields;

    bool m_colData_b = false;
    bool m_faces_b = false;
    bool m_mid_planes_b = false;
    bool m_project_coaption_b = false;
    bool m_project_main_direction_b = false;
    bool m_edges_b = false;
    bool m_boundary_edges_b = false;
    bool m_borders_b = false;
    bool m_leaf_profiles_b = false;

public:
    CoaptHeight(nets_t valve, InputProcessor gPrms, double margin);
    ~CoaptHeight();
    void computeCoaptationField();
    void clearCoaptationField();
    void savetoOFF(string directory);
    void savetoSTL(string directory, string name = "");
    void savetoSTL_remap(net_t& net, int id, string directory, string name = "");
    void saveCurrentMeshToSTL(string directory, string name = "mesh");
    void print_distribution(string directory, string name = "distribution", int N = 100);
    void print_billowing(string directory, string name = "billowing");
    double get_coapt_field_width(const pair<int, int>& id);
    vector<double> get_distribution(int N, const pair<int, int>& id);
    pair<vector<array<double, 2>>, double> get_distribution_per_leaf(int N, int mesh_n);
    void print_distribution_per_leaf(string directory, string name = "distribution_per_leaf", int N = 200);
    auto& get_billow_plate();
    double get_billowing_at(int nleaf);
    array<double, 3> get_billowing();
    nets_t get_mesh(){
        return mesh;
    }
    void refine_boundary_elems();
    void dump();
    void dump_planes();
    void test_to2d();
};

inline int testMeassure(InputProcessor& gPrms){
    nets_t valve = download_nets_from_file((gPrms.ro.res_dir + "result.nts").c_str());
    //to_stl(valve, (gPrms.ro.res_dir + "downloaded").c_str());

    CoaptHeight c(valve, gPrms, 0.11*3);
    //c.compute_plane_leaf_coaptation_profile(0);
    c.print_billowing(gPrms.ro.res_dir, "billowing");
    c.print_distribution(gPrms.ro.res_dir, "distribution-new", 100);
    c.refine_boundary_elems();
    c.refine_boundary_elems();
    c.refine_boundary_elems();
    c.print_distribution_per_leaf(gPrms.ro.res_dir, "distribution_per_leaf", 1000);
    c.savetoSTL(gPrms.ro.res_dir);
//    nets_t leaflet[3];
//    for (int i = 0; i < 3; ++i) {
//        leaflet[i] = formated_in(gPrms.l_ins[i].leaflet_file.c_str());
//        c.savetoSTL_remap(leaflet[i].nets[0], i, gPrms.ro.res_dir);
//    }
    //c.find_refined_elems();
    /*c.refine_boundary_elems();
    c.refine_boundary_elems();
    c.refine_boundary_elems();
    c.savetoOFF(gPrms.ro.res_dir);
    c.savetoSTL(gPrms.ro.res_dir);
    c.saveCurrentMeshToSTL(gPrms.ro.res_dir);
    c.print_distribution(gPrms.ro.res_dir, "distribution1", 1000);
    c.print_billowing(gPrms.ro.res_dir, "billowing1");*/

    return 0;
}

#endif