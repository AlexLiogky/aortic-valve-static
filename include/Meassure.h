#ifndef MEASSURE_H_INCLUDED
#define MEASSURE_H_INCLUDED

#include "nets.h"
#include "save-data.h"
#include "precomputation.h"
#include "World.h"
#include "InputProcessor.h"
#include <string>
#include <sstream>
#include <fstream>
#include <iostream>
#include <map>
#include <set>
#include <algorithm>
#include <cmath>

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
    void _savetoSTL(std::ofstream& stl, ElemIterator eit, ElemToPoints e_to_pnts, string name = ""){
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
    double field_directed_diam(const pair<int, int>& id, const Line_2& line, const Vector_2& dir);
    void compute_billow_plate();
    void find_refined_elems();
    vector<elem_t*> refine_elem(int mesh_n, elem_t* e);
    set<int> test_node_in_springs(node_t* n, net_t net);
    void compute_billowing_pnts();

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
    map<pair<int, int>, vector<Point_2>> projection;
    map<pair<int, int>, set<pair<int, int>>> m_edges;
    map<pair<int, int>, set<pair<int, int>>> m_boundary_edges;
    map<pair<int, int>, K::Vector_2> proj_dir;
    map<pair<int, int>, vector<array<int, 3>>> m_faces;
    map<pair<int, int>, pair<Point_2, Point_2>> borders;
    vector<vector<elem_t*>> m_no_coapt_elems;
    vector<map<elem_t*, bool>> m_refined_elems;
    RefineElem m_refiner;
    vector<map<int, pair<int, int>>> remap = {{{0, {0, 1}}, {1, {0, 2}}},
                                              {{0, {1, 0}}, {1, {1, 2}}},
                                              {{0, {2, 0}}, {1, {2, 1}}}};
    map<pair<int, int>, Plane> m_planes;
    To_colData _to_colData;
    point_t axis;

    bool m_colData_b = false;
    bool m_faces_b = false;
    bool m_mid_planes_b = false;
    bool m_project_coaption_b = false;
    bool m_project_main_direction_b = false;
    bool m_edges_b = false;
    bool m_boundary_edges_b = false;
    bool m_borders_b = false;

public:
    CoaptHeight(nets_t valve, InputProcessor gPrms, double margin);
    ~CoaptHeight();
    void computeCoaptationField();
    void clearCoaptationField();
    void savetoOFF(string directory);
    void savetoSTL(string directory, string name = "");
    void saveCurrentMeshToSTL(string directory, string name = "mesh");
    void print_distribution(string directory, string name = "distribution", int N = 100);
    void print_billowing(string directory, string name = "billowing");
    double get_coapt_field_width(const pair<int, int>& id);
    vector<double> get_distribution(int N, const pair<int, int>& id);
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
    to_stl(valve, (gPrms.ro.res_dir + "downloaded").c_str());

    CoaptHeight c(valve, gPrms, 0.11*2);
    c.print_billowing(gPrms.ro.res_dir, "billowing");
    c.print_distribution(gPrms.ro.res_dir, "distribution", 1000);
    //c.find_refined_elems();
    c.refine_boundary_elems();
    c.refine_boundary_elems();
    c.refine_boundary_elems();
    c.savetoOFF(gPrms.ro.res_dir);
    c.savetoSTL(gPrms.ro.res_dir);
    c.saveCurrentMeshToSTL(gPrms.ro.res_dir);
    c.print_distribution(gPrms.ro.res_dir, "distribution1", 1000);
    c.print_billowing(gPrms.ro.res_dir, "billowing1");

    return 0;
}

#endif