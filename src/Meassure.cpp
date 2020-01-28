#include <Meassure.h>

void CoaptHeight::set_colData(){
    if (m_colData_b) return;
    m_colData_b = true;
    nets_t empty_mesh = nets_t_get_net(0);
    World s(mesh, empty_mesh, _to_colData.conditions, _to_colData.solver_data, 0, 0, _to_colData.colission_margin);
    World::ColissionType data = s.getCollision(_to_colData.computive_mrg);

    map<pair<int, int>, set<node_t*>>& _colData = colData;
    for (auto& i: data){
        set<node_t*> v;
        v.insert(get<2>(i));
        _colData.insert(make_pair(make_pair(get<0>(i), get<1>(i)), v));
        _colData[make_pair(get<0>(i), get<1>(i))].insert(CoaptHeight::mesh.nets[get<0>(i)].vrtx.nodes[get<2>(i)->id]);
    }
    for (auto& i: _colData) {
        m_nodes.insert({i.first, vector<node_t*>()});
        auto& v = m_nodes[i.first];
        v.reserve(i.second.size());
        int counter = 0;
        for (auto &j: i.second) {
            if (!j->contact_elem_id){
                j->contact_elem_id = (int*)calloc(2, sizeof(int));
                j->contact_elem_id[0] = j->contact_elem_id[1] = -1;
            }
            j->contact_elem_id[mapre[i.first]] = counter++;
            v.push_back(j);
        }
    }
}


void CoaptHeight::computeCoaptationField(){
    set_colData();
    computeFaces();
    int state = 0;
    state = set_mid_planes();
    state = project_coaption_to_planes(state);
    project_main_direction(state);;
    set_all_edges();
    //set_triangle_edges();
    compute_boundary_edges();
    set_borders();
}

void CoaptHeight::clearCoaptationField(){
    colData.clear();
    m_nodes.clear();
    m_projection.clear();
    m_edges.clear();
    m_boundary_edges.clear();
    proj_dir.clear();
    m_faces.clear();
    m_borders.clear();
    m_no_coapt_elems.clear();
    m_refined_elems.clear();
    m_planes.clear();
    m_colData_b = m_faces_b = m_mid_planes_b = false;
    m_project_coaption_b = m_project_main_direction_b = false;
    m_edges_b = m_boundary_edges_b = m_borders_b = false;
    m_leaf_profiles_b = false;
}

CoaptHeight::CoaptHeight(nets_t valve, InputProcessor gPrms, double margin):
mesh{cpy_nets(valve)},
_to_colData {{gPrms.sp.delta, gPrms.sp.eps, gPrms.ots[0].elastic_model_type},
{gPrms.ots[0].pressure},
gPrms.ots[0].collision_margin,
margin
}
{
    auto& bl_flow = gPrms.a_in.main_flow_direction;
    axis = point_t_get_point(-bl_flow[0], -bl_flow[1], -bl_flow[2]);
}

CoaptHeight::~CoaptHeight(){
    nets_t_destruct(mesh);
}

void CoaptHeight::dump(){
    set_colData();
    for (auto& i: colData){
        std::cout << i.first.first << " -> " << i.first.second << ":\n";
        for (auto& j: i.second) {
            std::cout << "      ";
            point_t_dump(j->coord);
        }
    }
}

void CoaptHeight::savetoOFF(string directory){
    set_colData();
    for (auto& i: m_nodes){
        string name = directory + to_string(i.first.first) + "-" + to_string(i.first.second) + ".off";
        std::ofstream off_file(name);
        off_file << "OFF\n\n";
        off_file << i.second.size() <<" "<< m_faces[i.first].size() << " 0\n";
        for (auto& j: i.second)
        {
            off_file << j->coord.coord[0] << " " << j->coord.coord[1] << " " << j->coord.coord[2] << "\n";
        }
        for (auto& j: m_faces[i.first]){
            off_file << j.size() << "   " << j[0] << " " << j[1] << " " << j[2] << "\n";
        }
        off_file.close();
    }
}

void CoaptHeight::savetoSTL_remap(net_t& net, int id, string directory, string name){
    set_colData();
    computeFaces();

    set<int> saved_elems;
    for (int j = 0; j < mesh.nets[id].elems.count; ++j){
        elem_t& e = *mesh.nets[id].elems.elems[j];
        for (int k = 0; k < 2; ++k){
            set<node_t*>& nodes = colData[remap[id][k]];
            if (auto n0 = nodes.find(e.vrts[0]); n0 == end(colData[remap[id][k]])) continue;
            if (auto n1 = nodes.find(e.vrts[1]); n1 == end(colData[remap[id][k]])) continue;
            if (auto n2 = nodes.find(e.vrts[2]); n2 == end(colData[remap[id][k]])) continue;
            saved_elems.insert(e.id);
            }
    }
    string fname = directory + "coapt-" + to_string(id) + "-remap" + name + ".stl";
    std::ofstream stl(fname);
    auto e_to_pnts = [&net](auto& e){
        auto& n = net.elems.elems[e]->vrts;
        return array{n[0]->coord, n[1]->coord, n[2]->coord};
    };
    _savetoSTL(stl, saved_elems, e_to_pnts, name);
    stl.close();
    fname = directory + "nocoapt-" + to_string(id) + "-remap" + name + ".stl";
    stl.open(fname);
    auto ee_to_pnts = [&net](auto& e){
        auto& n = net.elems.elems[e->id]->vrts;
        return array{n[0]->coord, n[1]->coord, n[2]->coord};
    };
    _savetoSTL(stl, m_no_coapt_elems[id], ee_to_pnts, name);
}

void CoaptHeight::savetoSTL(string directory, string name){
    set_colData();
    computeFaces();
    for (auto& i: m_faces){
        string fname = directory + to_string(i.first.first) + "-" + to_string(i.first.second) + ".stl";
        std::ofstream stl(fname);
        auto e_to_pnts = [&](auto& e){
            auto& n = m_nodes[i.first];
            return array{n[e[0]]->coord, n[e[1]]->coord, n[e[2]]->coord};
        };
        _savetoSTL(stl, i.second, e_to_pnts);
    }
    for (int i = 0; i < m_no_coapt_elems.size(); ++i){
        string fname = directory + "no_coapt_" + to_string(i) + ".stl";
        std::ofstream stl(fname);
        auto e_to_pnts = [&](auto& e){
            auto n = e->vrts;
            return array{n[0]->coord, n[1]->coord, n[2]->coord};
        };
        _savetoSTL(stl, m_no_coapt_elems[i], e_to_pnts);
    }
    set_colData();
    computeFaces();
    int state = 0;
    state = set_mid_planes();
    state = project_coaption_to_planes(state);
    for (auto& i: m_faces){
        string fname = directory + to_string(i.first.first) + "-" + to_string(i.first.second) + "_projected" +".stl";
        std::ofstream stl(fname);
        //cout << "saveSTL: p: " << m_planes[i.first].point() << endl;
        auto P2p = [&](const Point_2& p) {
            auto& pl = m_planes[i.first];
            auto vv = pl.point() + (p[0] / sqrt(pl.base1().squared_length()) * pl.base1() + p[1] / sqrt(pl.base2().squared_length()) * pl.base2());
            return VEC(vv[0], vv[1], vv[2]);
        };
        auto e_to_pnts = [&](auto& e){
            auto& n = m_projection[i.first];
            return array{P2p(n[e[0]]), P2p(n[e[1]]), P2p(n[e[2]])};
        };
        _savetoSTL(stl, i.second, e_to_pnts);
    }

    compute_leaf_profiles();
    auto nconv1 = [](node_t* n) {return n->id;};
    auto invmap = [](map<int, int> m){
        map<int, int> res;
        for (const auto& i: m)
            res.insert({i.second, i.first});
        return res;
    };
    for (auto it = m_leaflet_coapt_fields.begin(); it != m_leaflet_coapt_fields.end(); ++it) {
        string fname = directory + "coapt_profile_" + to_string(it->first + 1) + ".stl";
        std::ofstream stl(fname);
        auto& [pnts, links, gmap, borders] = it->second;
        auto check_link = [&links = links, inmp = invmap(gmap)](pair<int, int> q) mutable {
            int tmp[2] = {inmp[q.first], inmp[q.second]};
            if (tmp[0] > tmp[1]) swap(tmp[0], tmp[1]);
            return (links.count({tmp[0], tmp[1]}) > 0);
        };
        auto inmp = invmap(gmap);
        auto e_to_pnts3 = [&pnts = pnts, &inmp = inmp](auto &e) {
            array<point_t, 3> arr;
            for (int i = 0; i < 3; ++i)
                arr[i] = VEC(pnts[inmp[e[i]]][0], pnts[inmp[e[i]]][1], 0);
            return arr;
        };
        auto triags = convert_edges_to_elems(nconv1, check_link, it->first);
        _savetoSTL(stl, triags, e_to_pnts3);
    }
}

void CoaptHeight::computeFaces(){
    if (m_faces_b) return;
    m_faces_b = true;
    m_no_coapt_elems.resize(mesh.count);
    for (int i = 0; i < mesh.count; ++i){
        m_no_coapt_elems[i].reserve(mesh.nets[i].elems.count);
        for (int j = 0; j < mesh.nets[i].elems.count; ++j){
            elem_t& e = *mesh.nets[i].elems.elems[j];
            bool not_coapt = true;
            for (int k = 0; k < 2; ++k){
                set<node_t*>& nodes = colData[remap[i][k]];
                m_faces.insert({remap[i][k], vector<array<int, 3>>()});
                auto n0 = nodes.find(e.vrts[0]);
                if (n0 == end(colData[remap[i][k]])) continue;
                auto n1 = nodes.find(e.vrts[1]);
                if (n1 == end(colData[remap[i][k]])) continue;
                auto n2 = nodes.find(e.vrts[2]);
                if (n2 == end(colData[remap[i][k]])) continue;
                not_coapt = false;
                if (e.coef > 0)
                    m_faces[remap[i][k]].push_back({(int)(*n1)->contact_elem_id[k], (int)(*n0)->contact_elem_id[k], (int)(*n2)->contact_elem_id[k]});
                else
                    m_faces[remap[i][k]].push_back({(int)(*n0)->contact_elem_id[k], (int)(*n1)->contact_elem_id[k], (int)(*n2)->contact_elem_id[k]});
            }
            if (not_coapt)
                m_no_coapt_elems[i].push_back(&e);
        }
    }
}

CoaptHeight::Point_3 CoaptHeight::point_t_to_Point(const point_t& p0) const {
    auto* p = p0.coord;
    return Point_3(p[0], p[1], p[2]);
}

map<pair<int, int>, CoaptHeight::Plane> CoaptHeight::find_mid_planes(const int dim) {
    map<pair<int, int>, Plane> planes;
    auto to_pnt = [](node_t *n) {
        double *p = n->coord.coord;
        return Point_3(p[0], p[1], p[2]);
    };
    auto convT = [=](array<node_t *, 3> tr) {
        return Triangle(to_pnt(tr[0]), to_pnt(tr[1]), to_pnt(tr[2]));
    };
    switch (dim) {
        case 0: {
            auto &refs = m_nodes;
            for (int i = 0, cnt = refs.size() / 2; i < cnt; ++i) {
                for (int j = i + 1; j < cnt; ++j) {
                    Plane plane = get_linear_least_squares_fitting_plane<0>(refs[{i, j}],refs[{j, i}],
                                                                            to_pnt);
                    planes.insert({{i, j}, plane});
                    planes.insert({{j, i}, plane});
                }
            }
        }
        case 1 ... 2: {
            map<pair<int, int>, vector<array<node_t *, 3>>> refs;
            for (auto &i: m_faces) {
                refs.insert({i.first, vector<array<node_t *, 3>>()});
                auto &v = refs[i.first];
                v.reserve(i.second.size());
                auto &nd = m_nodes[i.first];
                for (auto &j: i.second) {
                    v.push_back({nd[j[0]], nd[j[1]], nd[j[2]]});
                }
            }

            for (int i = 0, cnt = refs.size() / 2; i < cnt; ++i)
                for (int j = i + 1; j < cnt; ++j) {
                    Plane plane;
                    if (dim == 2)
                        plane = get_linear_least_squares_fitting_plane<2>(refs[{i, j}],refs[{j, i}], convT);
                    else
                        plane = get_linear_least_squares_fitting_plane<1>(refs[{i, j}],refs[{j, i}], convT);
                    planes.insert({{i, j}, plane});
                    planes.insert({{j, i}, plane});
                }

            return planes;
        }

        default: return planes;
    }

}

int CoaptHeight::set_mid_planes(const int dim){
    static int loc_dim = -1;
    if (m_mid_planes_b && loc_dim == dim)
        return 1;
    m_mid_planes_b = true;
    loc_dim = dim;
    m_planes = find_mid_planes(dim);
    if (mesh.count == 3){
        double dd[3][2] = {{1, 2}, {0, 2}, {0, 1}};
        int changed[3] = {0, 0, 0};
        for (int i = 0; i < 3; ++i){
            Plane   &p1 = m_planes[remap[i][0]],
                    &p2 = m_planes[remap[i][1]];
            if (CGAL::scalar_product(p1.orthogonal_vector(), p2.orthogonal_vector()) < 0){
                int to = remap[i][0].first + remap[i][0].second - 1;
                if (changed[to]){
                    Plane &pp = p2;
                    Vector_3 v = pp.orthogonal_vector();
                    pp = Plane(-pp.a(), -pp.b(), -pp.c(), -pp.d());
                    to = remap[i][1].first + remap[i][1].second - 1;
                    changed[to]++;
                }
                else{
                    Plane &pp = p1;
                    Vector_3 v = pp.orthogonal_vector();
                    pp = Plane(-pp.a(), -pp.b(), -pp.c(), -pp.d());
                    changed[to]++;
                }
            }
        }
    }
    return 0;
}

void CoaptHeight::dump_planes(){
    set_colData();
    computeFaces();
    int state = 0;
    state = set_mid_planes();
    for (auto& i: m_planes)
        cout << i.first.first << " - " << i.first.second << " : " << i.second << endl;
}

void CoaptHeight::test_to2d(){
    set_colData();
    computeFaces();
    int state = 0;
    state = set_mid_planes();
    pair id = {0, 1};
    for (int k = 0; k < 5; ++k) {
        point_t p = m_nodes[id][k]->coord;
        auto plane = m_planes[id];
        Point_3 p_3 = plane.projection(point_t_to_Point(p));
        Vector_3 v_3 = plane.projection(point_t_to_Point(p)) - plane.point();
        Point_2 p_2 = plane.to_2d(point_t_to_Point(p));
        Point_2 _p_2(CGAL::scalar_product(v_3, plane.base1()), CGAL::scalar_product(v_3, plane.base2()));
        cout << "base1: " << plane.base1() << " base2: " << plane.base2() << endl;
        cout << "to_2d: " << p_2 << endl;
        cout << "projection to base: " << _p_2 << endl;
        cout << "diff = " << _p_2 - p_2 << endl << endl;
    }
}

CoaptHeight::Point_2 CoaptHeight::p2plane(const point_t& p, const Plane& plane){
    auto pp = plane.point();
    array v = {plane.base1() / sqrt(plane.base1().squared_length()), plane.base2() / sqrt(plane.base2().squared_length())};
    auto ref_p = point_t_to_Point(p) - pp;
    return {CGAL::scalar_product(ref_p, v[0]), CGAL::scalar_product(ref_p, v[1])};
}

int CoaptHeight::project_coaption_to_planes(int state){
    if (state && m_project_coaption_b) return 1;
    m_project_coaption_b = true;
    for (auto& i: m_nodes){
        m_projection.insert({i.first, vector<Point_2>()});
        auto& pnts = m_projection[i.first];
        pnts.reserve(i.second.size());
        auto& plane = m_planes[i.first];
        auto p = plane.point();
        array v = {plane.base1(), plane.base2()};
        //cout << "project_coaption_to_planes: p: " << plane.point() << endl;
        for (auto& j: i.second){
            auto ref_p = point_t_to_Point(j->coord) - p;
            pnts.push_back({CGAL::scalar_product(ref_p, v[0])/sqrt(v[0].squared_length()), CGAL::scalar_product(ref_p, v[1])/sqrt(v[1].squared_length())});
        }
    }
    return 0;
}

int CoaptHeight::project_main_direction(int state){
    if (state && m_project_main_direction_b) return 1;
    m_project_main_direction_b = true;
    auto& direction = axis;
    K::Vector_3 v{direction.coord[0], direction.coord[1], direction.coord[2]};
    for (auto& i: m_planes){
        K::Vector_2 dir = {CGAL::scalar_product(v, i.second.base1()), CGAL::scalar_product(v, i.second.base2())};
        dir /= sqrt(dir.squared_length());
        proj_dir.insert({i.first, dir});
    }
    return 0;
}


pair<CoaptHeight::Point_2, CoaptHeight::Point_2> CoaptHeight::get_minmax_Point_set(vector<CoaptHeight::Point_2>& data, CoaptHeight::Vector_2 dir){
    auto comp = [=](const Point_2& p1, const Point_2& p2){
        Vector_2 v1(p1[0], p1[1]);
        Vector_2 v2(p2[0], p2[1]);
        return CGAL::scalar_product(dir, v1) < CGAL::scalar_product(dir, v2);
    };
    auto mm = minmax_element(data.begin(), data.end(), comp);
    return {*(mm.first), *(mm.second)};
}


void CoaptHeight::set_borders(){
    if (m_borders_b) return;
    m_borders_b = true;
    for (auto& i: m_projection){
        Vector_2 dir(proj_dir[i.first][1], -proj_dir[i.first][0]);
        m_borders.insert({i.first, get_minmax_Point_set(i.second, dir)});
    }

    /*auto pvt = [](Point_2 p) { return Vector_2(p[0], p[1]); };
    for (int i = 0; i < m_nodes.size(); ++i)
        for (int j = i + 1; j < m_nodes.size(); ++j){
            Vector_2 dir(proj_dir[{i, j}][1], -proj_dir[{i, j}][0]);
            auto p0 = m_borders[{i, j}].first, p1 = m_borders[{i, j}].second;
            if (dir * pvt(p0) > dir * pvt(m_borders[{j, i}].first)) p0 = m_borders[{j, i}].first;
            if (dir * pvt(p1) < dir * pvt(m_borders[{j, i}].second)) p1 = m_borders[{j, i}].second;
            m_borders[{i, j}] = m_borders[{j, i}] = {p0, p1};
        }*/
}


void CoaptHeight::print_distribution(string directory, string name, int N){
    computeCoaptationField();
    string fname = directory + name + ".csv";
    ofstream csv(fname);
    pair<int, int> right_map[] = {{0,1}, {1,0}, {0, 2}, {2, 0}, {1, 2}, {2, 1}};
    for (int j = 0; j < m_boundary_edges.size(); ++j){
        vector<double> distr = get_distribution(N, right_map[j]);
        double width = get_coapt_field_width(right_map[j]);
        csv << right_map[j].first + 1 << " -> " << right_map[j].second + 1 << "; width = " << width << ", ";
//            cout << "borders: " << borders[i.first].first << ", " << borders[i.first].second << endl;
        for (auto& i: distr){
            csv << i << ", ";
        }
        csv << endl;
    }
}

void CoaptHeight::set_triangle_edges(){
    for (auto& i: m_faces){
        m_edges.insert({i.first, set<pair<int, int>>()});
        auto& edges = m_edges[i.first];
        for (auto j: i.second){
            sort(j.begin(), j.end());
            edges.insert({j[0], j[1]});
            edges.insert({j[0], j[2]});
            edges.insert({j[1], j[2]});
        }
        //cout << "3 * trs = " << 3 * i.second.size() << " edges = " << edges.size() << endl;
    }
}

void CoaptHeight::set_all_edges() {
    if (m_edges_b) return;
    m_edges_b = true;
    for (int i = 0; i < mesh.count; ++i) {
        for (int j = 0; j < mesh.nets[i].springs.count; ++j) {
            spring_t &e = *mesh.nets[i].springs.springs[j];
            for (int k = 0; k < 2; ++k) {
                set<node_t *> &nodes = colData[remap[i][k]];
                m_edges.insert({remap[i][k], set<pair<int, int>>()});
                auto n0 = nodes.find(e.ends[0]);
                if (n0 == end(colData[remap[i][k]])) continue;
                auto n1 = nodes.find(e.ends[1]);
                if (n1 == end(colData[remap[i][k]])) continue;
                int id[2] = {(int) (*n0)->contact_elem_id[k], (int) (*n1)->contact_elem_id[k]};
                sort(id, id + 2);
                m_edges[remap[i][k]].insert({id[0], id[1]});
            }
        }
    }
}

void CoaptHeight::compute_boundary_edges(){
    if (m_boundary_edges_b) return;
    m_boundary_edges_b = true;
    for (auto& i: m_edges){
        map<pair<int, int>, int> in_elems;
        for (auto j: m_faces[i.first]){
            sort(j.begin(), j.end());
            for (int ii = 0; ii < 3; ++ii)
                for (int jj = ii + 1; jj < 3; ++jj) {
                    in_elems.insert({{j[ii], j[jj]}, 0});
                    ++in_elems[{j[ii], j[jj]}];
                }
        }
        m_boundary_edges.insert({i.first, set<pair<int, int>>()});
        auto& data = m_boundary_edges[i.first];
        for (auto& j: i.second){
            auto it = in_elems.find(j);
            if (it == in_elems.end() || it->second < 2)
                data.insert(j);
        }
    }
}

array<double, 2> CoaptHeight::_clever_field_directed_diam(const vector<Point_2>& pnts, const set<pair<int, int>>& edges, const Line_2& line, const Vector_2& dir){
    double mmin = 1e20, mmax = -1e20;
    for (auto& i: edges){
        Segment_2  edge(pnts[i.first], pnts[i.second]);
        auto result = CGAL::intersection(line, edge);
        Point_2 p;
        if (result) {
            if (const Segment_2* s = boost::get<Segment_2>(&*result)) {
                p = (*s).source();
            } else {
                const Point_2* _p = boost::get<Point_2 >(&*result);
                p = *_p;
            }
        }
        else continue;
        Vector_2 pp(p[0], p[1]);
        double res = pp * dir;
        mmin = min(mmin, res);
        mmax = max(mmax, res);
    }
    return {mmin, mmax};
}

double CoaptHeight::_field_directed_diam(const vector<Point_2>& pnts, const set<pair<int, int>>& edges, const Line_2& line, const Vector_2& dir){
    auto diam = _clever_field_directed_diam(pnts, edges, line, dir);
    double distance = diam[1] - diam[0];
    return (distance > 0) ? distance : 0;
}

double CoaptHeight::field_directed_diam(const pair<int, int>& id, const Line_2& line, const Vector_2& dir){
    return _field_directed_diam(m_projection[id],  m_boundary_edges[id], line, dir);
}

double CoaptHeight::get_coapt_field_width(const pair<int, int>& id){
    computeCoaptationField();
    auto _brds = m_borders[id];
    Vector_2 dir = proj_dir[id];
    Vector_2 dir_p(dir[1], -dir[0]);
    Vector_2 p[2] = {{_brds.first[0], _brds.first[1]}, {_brds.second[0], _brds.second[1]}};
    return  fabs((p[1] - p[0])*dir_p);
}
vector<double> CoaptHeight::get_distribution(int N, const pair<int, int>& id){
    computeCoaptationField();
    Point_2 _p_min =  get_minmax_Point_set(m_projection[id], proj_dir[id]).first;
    Vector_2 p_min(_p_min[0], _p_min[1]);
    Vector_2 dir = proj_dir[id];
    Vector_2 dir_p(dir[1], -dir[0]);
    auto _brds = m_borders[id];
    Vector_2 p[2] = {{_brds.first[0], _brds.first[1]}, {_brds.second[0], _brds.second[1]}};
    vector<double> distribution;
    distribution.reserve(N+1);
    for (int i = 0; i < N+1; ++i){
        double cf = static_cast<double>(i) / N;
        const double eps = 1e-5;
        auto point = (p_min * dir - eps) * dir + ((p[0] + (cf * (p[1] - p[0])))*dir_p)*dir_p;
        Point_2 l_p(point[0], point[1]);
        distribution.push_back(field_directed_diam(id, Line_2(l_p, dir), dir));
    }
    return distribution;
}

void CoaptHeight::compute_billow_plate(){
    for (int i = 0; i < mesh.count; ++i){
        auto& nodes = mesh.nets[i].vrtx;
        point_t p;
        double last_mul = 0;
        bool finded = false;
        for (int j = 0; j < nodes.count; ++j){
            if (is_fix(nodes.nodes[j]->state)){
                if (finded){
                    double mul = DOT(nodes.nodes[j]->coord, axis);
                    if (mul < last_mul){
                        last_mul = mul;
                        p = nodes.nodes[j]->coord;
                    }
                }
                else {
                    p = nodes.nodes[j]->coord;
                    last_mul = DOT(p, axis);
                    finded = true;
                }
            }
        }
        m_billow_plate[i] = p;
    }
    m_billow_plate_state = true;
}

auto& CoaptHeight::get_billow_plate(){
    if (m_billow_plate_state) return m_billow_plate;
    compute_billow_plate();
    return m_billow_plate;
}

void CoaptHeight::find_refined_elems(){
    set_colData();
    computeFaces();
    vector<map<elem_t*, bool>>& refined = m_refined_elems;
    refined.resize(mesh.count);
    for (int i = 0; i < mesh.count; ++i)
        for (int j = 0; j < mesh.count; ++j){
            if (j == i) continue;
            set<node_t*> data;
            auto& init = m_nodes[{i, j}];
            data.insert(init.begin(), init.end());
            for (auto& k: m_no_coapt_elems[i]){
                if (data.count(k->vrts[0]) || data.count(k->vrts[1]) || data.count(k->vrts[2]))
                    refined[i].insert({k, false});
            }
        }
}


    auto CoaptHeight::RefineElem::get_other_node(elem_t* elem, node_t* n1, node_t* n2){
        if (!elem) return (node_t*)nullptr;
        node_t* next = elem->vrts[0];
        if (next == n1 || next == n2) next = elem->vrts[1];
        if (next == n1 || next == n2) next = elem->vrts[2];
        return next;
    }

    auto CoaptHeight::RefineElem::get_neigh_elem(const net_t net, const elem_t* elem, const node_t* n1, const node_t* n2){
        const node_t* from[2] = {n1, n2};
        elem_t* res = nullptr;
        for (int j = 0; j < elem->cnt_neighbours; ++j){
            set<const node_t*> nds;
            elem_t& el = *net.elems.elems[elem->neighbours_id[j]];
            nds.insert(el.vrts, el.vrts + 3);
            if (nds.count(from[0]) && nds.count(from[1]) && &el != elem){
                res = &el;
                continue;
            }
        }
        return res;
    }

    pair<array<elem_t*, 3>, int> CoaptHeight::RefineElem::compute_shared_elems(net_t net, elem_t* e){
        array<elem_t*, 3> shared_elems;
        int shared_elems_cnt = 0;
        for (int i = 0; i < 3; ++i){
            auto elem = get_neigh_elem(net, e, e->vrts[i %3], e->vrts[(i + 1) % 3]);
            shared_elems[i] = elem;
            if (elem){
                ++shared_elems_cnt;
            }
        }
        return {shared_elems, shared_elems_cnt};
    }

    int CoaptHeight::RefineElem::get_state(int state1, int state2) {
        if (state1 == state2) return state1;
        if (is_fix(state1) && is_fix(state2))
            return FIX_BND;
        return IN;
    }


    elem_t* CoaptHeight::RefineElem::e_construct(node_t* node1, node_t* node2, node_t* node3, unsigned int id, elem_t* on){
        m_nds_els_id[node1->id].insert(id);
        m_nds_els_id[node2->id].insert(id);
        m_nds_els_id[node3->id].insert(id);
        elem_t* elem = elem_t_construct(node1, node2, node3, id);
        elem->coef *= (DOT(elem->or_area, on->or_area) < 0) ? -1 : 1;
        return elem;
    }
    spring_t* CoaptHeight::RefineElem::spr_construct(node_t* node1, node_t* node2, unsigned int id){
        m_nds_sprs_id[node1->id].insert(id);
        m_nds_sprs_id[node2->id].insert(id);
        return spring_t_construct(node1, node2, id);
    }

    void CoaptHeight::RefineElem::set_data_to_existing_node_0(node_t* n, set<int>& excluding, set<int> what)
    {
        auto& e_ids = excluding;
        auto& data = what;
        for (int j = 0; j < n->cnt_elems; ++j) {
            int loc_id = n->elems_id[j];
            if (e_ids.count(loc_id) == 0 /*&& loc_id < m_net.elems.count*/)
                data.insert(loc_id);
        }
        n->cnt_elems   = data.size();
        n->elems_id    = (unsigned int*) realloc(n->elems_id, n->cnt_elems * sizeof(unsigned int));
        copy(data.begin(), data.end(), n->elems_id);
    }
    void CoaptHeight::RefineElem::set_data_to_existing_node_1(node_t* n, set<int>& excluding, set<int> what){
        auto &s_ids = excluding;
        auto &sdata = what;
        for (int j = 0; j < n->cnt_springs; ++j) {
            int loc_id = n->springs_id[j];
            if (s_ids.count(loc_id) == 0 /*&& loc_id < m_net.springs.count*/)
                sdata.insert(loc_id);
        }
        n->cnt_springs = sdata.size();
        n->springs_id = (unsigned int *) realloc(n->springs_id, n->cnt_springs * sizeof(unsigned int));
        copy(sdata.begin(), sdata.end(), n->springs_id);
    }

    void CoaptHeight::RefineElem::compute_new_indeces(){
        auto data = compute_shared_elems(*m_net, m_e);
        for (int i = 0; i < 3; ++i)
            if (data.first[i])
                m_neigh_elems.insert({i, {data.first[i], get_other_node(data.first[i], m_e->vrts[i % 3], m_e->vrts[(i+1)%3])}});
        int shr_cnt = m_neigh_elems.size();
        m_nodes_id.resize(6 + shr_cnt);
        m_springs_id.resize(9 + shr_cnt);
        m_elems_id.resize(4 + 2 * shr_cnt);
        m_elems_id[0] = m_e->id;
        for (int i = 0; i < 3; ++i){
            m_nodes_id[i] = (int)m_e->vrts[i]->id;
            m_nodes_id[i + 3] = m_net->vrtx.count + i;
            m_springs_id[i] = get_shared_spring(*m_net, m_e->vrts[i % 3], m_e->vrts[(i + 1) % 3])->id;
            m_springs_id[i + 3] = m_net->springs.count + i;
            m_springs_id[i + 6] = m_net->springs.count + 3 + i;
            m_elems_id[1 + i] = m_net->elems.count + i;
        }
        for (int i = 0, j = 0; i < 3; ++i){
            node_t* node = get_other_node(get_neigh_elem(*m_net, m_e, m_e->vrts[i % 3], m_e->vrts[(i+1)%3]), m_e->vrts[i % 3], m_e->vrts[(i+1)%3]);
            if (node) {
                m_springs_id[j + 9] = m_net->springs.count + 6 + j;
                m_nodes_id[j + 6] = node->id;
                m_elems_id[4 + 2 * j] = data.first[i]->id;
                m_elems_id[4 + 2 * j + 1] = m_net->elems.count + 3 + j;
                ++j;
            }
        }
        for (auto& i: m_nodes_id){
            m_nds_els_id.insert({i, set<int>()});
            m_nds_sprs_id.insert({i, set<int>()});
        }
    }

    void CoaptHeight::RefineElem::create_new_objects(){
        m_new_nodes.reserve(3);
        for (int i = 0; i < 3; ++i){
            node_t* from[2] = {m_e->vrts[i % 3], m_e->vrts[(i + 1) % 3]};
            point_t pnt = SCAL_SUM(0.5, from[0]->coord, 0.5, from[1]->coord);
            double h = (from[0]->h + from[1]->h) / 2;
            int state = get_state(from[0]->state, from[1]->state);
            m_new_nodes.push_back(node_t_construct(pnt.coord[0], pnt.coord[1], pnt.coord[2], h, state, m_nodes_id[i + 3]));
            m_new_nodes[i]->initial = SCAL_SUM(0.5, from[0]->initial, 0.5, from[1]->initial);
        }

        m_new_elems.resize(1 + 3 + 2 * m_neigh_elems.size());
        m_new_elems[0] = e_construct(m_new_nodes[0], m_new_nodes[1], m_new_nodes[2], m_elems_id[0], m_e);
        for (int i = 0; i < 3; ++i){
            m_new_elems[i+1] = e_construct(m_e->vrts[i % 3], m_new_nodes[i % 3], m_new_nodes[(i + 2)%3], m_elems_id[i + 1], m_e);
        }
        int off = 4;
        for (auto& i: m_neigh_elems){
            m_new_elems[off] = e_construct(m_e->vrts[i.first % 3], m_new_nodes[i.first % 3], i.second.second, m_elems_id[off], i.second.first);
            off++;
            m_new_elems[off] = e_construct(m_e->vrts[(i.first + 1)% 3], m_new_nodes[i.first % 3], i.second.second, m_elems_id[off], i.second.first);
            off++;
        }

        m_new_sprs.resize(3 + 2 * 3 + m_neigh_elems.size());
        for (int i = 0; i < 3; ++i){
            m_new_sprs[i] = spr_construct(m_new_nodes[i % 3], m_new_nodes[(i + 1) % 3], m_springs_id[i]);
        }
        for (int i = 0; i < 3; ++i){
            m_new_sprs[2 * i + 3] = spr_construct(m_e->vrts[i % 3], m_new_nodes[i % 3], m_springs_id[3 + 2 * i]);
            m_new_sprs[2 * i + 1 + 3] = spr_construct(m_e->vrts[(i + 1) % 3], m_new_nodes[i % 3], m_springs_id[3 + 2 * i + 1]);
        }
        off = 9;
        for (auto& i: m_neigh_elems){
            m_new_sprs[off] = spr_construct(i.second.second, m_new_nodes[i.first], m_springs_id[off]);
            off++;
        }
    }
    void CoaptHeight::RefineElem::update_old_nodes(){
        for (int i = 0; i < 3; ++i){
            node_t* n = m_net->vrtx.nodes[m_nodes_id[i]];
            set<int> e_ids;
            e_ids.insert((int)m_e->id);
            if (auto it = m_neigh_elems.find(i); it != m_neigh_elems.end())
                e_ids.insert(it->second.first->id);
            if (auto it = m_neigh_elems.find((i + 2) % 3); it != m_neigh_elems.end())
                e_ids.insert(it->second.first->id);
            set_data_to_existing_node<0>(n, e_ids, m_nds_els_id[n->id]);
            set<int> s_ids;
            s_ids.insert(m_springs_id[i]);
            s_ids.insert(m_springs_id[(i + 2)%3]);
            set_data_to_existing_node<1>(n, s_ids, m_nds_sprs_id[n->id]);
        }
        for (auto& i: m_neigh_elems){
            auto& n = i.second.second;
            set<int> e_ids;
            e_ids.insert(i.second.first->id);
            set_data_to_existing_node<0>(n, e_ids, m_nds_els_id[n->id]);
            set<int> s_ids;
            set_data_to_existing_node<1>(n, s_ids, m_nds_sprs_id[n->id]);
        }
    }
    void CoaptHeight::RefineElem::insert_new_objs_into_net(){
        int off = 0;
        off = 3;
        m_net->vrtx.count += off;
        //m_net.vrtx.nodes = (node_t**)realloc(m_net.vrtx.nodes, m_net.vrtx.count * sizeof(node_t*));
        for (int i = 0; i < 3; ++i){
            set<int> e_ids;
            set_data_to_existing_node<0>(m_new_nodes[i], e_ids, m_nds_els_id[m_new_nodes[i]->id]);
            set_data_to_existing_node<1>(m_new_nodes[i], e_ids, m_nds_sprs_id[m_new_nodes[i]->id]);
            m_net->vrtx.nodes[m_new_nodes[i]->id] = m_new_nodes[i];
        }
        off = 3 + m_neigh_elems.size();
        m_net->elems.count += off;
        //m_net.elems.elems = (elem_t**)realloc(m_net.elems.elems, m_net.elems.count * sizeof(elem_t*));
        for (auto& i: m_new_elems) {
            if (i->id < m_net->elems.count - off)
                elem_t_destruct(m_net->elems.elems[i->id]);
            m_net->elems.elems[i->id] = i;
        }
        off = 6 + m_neigh_elems.size();
        m_net->springs.count += off;
        //m_net.springs.springs = (spring_t**)realloc(m_net.springs.springs, m_net.springs.count * sizeof(spring_t*));
        for (auto i: m_new_sprs){
            if (i->id < m_net->springs.count - off)
                spring_t_destruct(m_net->springs.springs[i->id]);
            m_net->springs.springs[i->id] = i;
        }

    }
    void CoaptHeight::RefineElem::prepare_elems(){
        set<int> changed_elems;
        for (auto& i: m_nodes_id){
            node_t* n = m_net->vrtx.nodes[i];
            for (int j = 0; j < n->cnt_elems; ++j) {
                changed_elems.insert(n->elems_id[j]);
            }
        }
        for (auto i: changed_elems){
            elem_t* e = m_net->elems.elems[i];
            if (e->neighbours_id){
                free(e->neighbours_id);
                e->cnt_neighbours = 0;
            }
            net_t_elem_t_set_elems_neighbours(*m_net, e);
        }
    }

    void CoaptHeight::RefineElem::clear(){
        if (!filled) return;
        m_net = nullptr;
        m_e = nullptr;
        m_neigh_elems.clear();
        m_nodes_id.resize(0);
        m_springs_id.resize(0);
        m_elems_id.resize(0);
        m_new_nodes.resize(0);
        m_new_sprs.resize(0);
        m_new_elems.resize(0);
        m_nds_els_id.clear();
        m_nds_sprs_id.clear();
    }

    void CoaptHeight::RefineElem::test(){
        for (int ii = 0; ii <  m_elems_id.size(); ++ii){
            int i = m_elems_id[ii];
            for (int j = 0; j < 3; ++j){
                elem_t* l_e = m_net->elems.elems[i];
                if (!get_shared_spring(*m_net, l_e->vrts[j % 3], l_e->vrts[(j + 1) % 3]))
                    cout << "Something break" << endl;
            }
        }
    }

    vector<elem_t*> CoaptHeight::RefineElem::refine_elem(net_t& net, elem_t* e){
        clear();
        m_net = &net; m_e = e; filled = true;
        compute_new_indeces();
        create_new_objects();
        update_old_nodes();
        insert_new_objs_into_net();
        prepare_elems();
        //test();
        vector<elem_t*> res;
        res.reserve(1 + m_neigh_elems.size());
        res.push_back(m_e);
        for(auto& i: m_neigh_elems)
            res.push_back(i.second.first);
        return res;
    }


vector<elem_t*> CoaptHeight::refine_elem(int mesh_n, elem_t* e){
    return m_refiner.refine_elem(mesh.nets[mesh_n], e);
}

set<int> CoaptHeight::test_node_in_springs(node_t* n, net_t net){
    set<int> res;
    for (int i = 0; i < net.springs.count; ++i)
        if (spring_t_node_belong(net.springs.springs[i], n))
            res.insert(i);
    return res;
}


void CoaptHeight::refine_boundary_elems(){
    find_refined_elems();
    for (int i = 0; i < m_refined_elems.size(); ++i){
        mesh.nets[i].vrtx.nodes = (node_t**) realloc(mesh.nets[i].vrtx.nodes, (mesh.nets[i].vrtx.count + mesh.nets[i].springs.count) * sizeof(node_t*));
        mesh.nets[i].springs.springs = (spring_t**) realloc(mesh.nets[i].springs.springs, 3 * mesh.nets[i].springs.count * sizeof(spring_t*));
        mesh.nets[i].elems.elems = (elem_t**) realloc(mesh.nets[i].elems.elems, 4 * mesh.nets[i].elems.count * sizeof(elem_t*));
        int off = 0;
        auto it = m_refined_elems[i].rbegin();
        for (int jj = 0; jj < m_refined_elems[i].size(); ++jj, ++it){
            auto& j = *it;
            if (!j.second){
                off++;
                auto res = refine_elem(i, j.first);
                for (auto& k: res)
                    if (auto it = m_refined_elems[i].find(k); it != m_refined_elems[i].end())
                        it->second = true;
            }
        }
    }
    clearCoaptationField();
}

void CoaptHeight::compute_billowing_pnts(){
    auto& plate = get_billow_plate();
    point_t normal = NORM(OR_AREA(plate[0], plate[1], plate[2]));
    if (DOT(normal, axis) < 0) SCAL_S(-1, &normal);
    for (int i = 0; i < mesh.count; ++i){
        auto& nodes = mesh.nets[i].vrtx;
        point_t p;
        double last_mul = 0;
        bool finded = false;
        for (int j = 0; j < nodes.count; ++j){
            if (finded){
                double mul = DOT(nodes.nodes[j]->coord, normal);
                if (mul < last_mul){
                    last_mul = mul;
                    p = nodes.nodes[j]->coord;
                }
            }
            else {
                p = nodes.nodes[j]->coord;
                last_mul = DOT(p, normal);
                finded = true;
            }
        }
        point_t& billow_pnt = p;
        m_billow_pnt[i] = billow_pnt;
        double distance = DOT(normal, DIF(billow_pnt, plate[0]));
        m_billowing[i] = distance;
    }
    m_billowing_state = true;
}

double CoaptHeight::get_billowing_at(int nleaf){
    if (nleaf >= 0 && nleaf < 3)
        return get_billowing()[nleaf];
    else return 1.0e20;
}

array<double, 3> CoaptHeight::get_billowing(){
    if (m_billowing_state) return m_billowing;
    compute_billowing_pnts();
    return m_billowing;
}

void CoaptHeight::saveCurrentMeshToSTL(string directory, string name){
    string fname = directory + name;
    to_stl(mesh, fname.c_str());
}

void CoaptHeight::print_billowing(string directory, string name){
    string fname = directory + name + ".txt";
    ofstream fout(fname);
    auto pnt2str = [](point_t& p){
        return "(" + to_string(p.coord[0]) + ", " + to_string(p.coord[1]) + ", " + to_string(p.coord[2]) + ")";
    };
    auto default_to_str = [](auto item){
        return to_string(item);
    };
    auto printer = [&](auto& container, string text, auto to_string){
        fout << text;
        for (auto& i: container){
            fout << to_string(i) << ", ";
        }
        fout << endl;
    };
    get_billowing();
    printer(m_billow_plate, "billow plate: ", pnt2str);
    printer(m_billow_pnt, "billow pnts: ", pnt2str);
    printer(m_billowing, "billowing: ", default_to_str);
}


void CoaptHeight::print_distribution_per_leaf(string directory, string name, int N){
    computeCoaptationField();
    string fname = directory + name + ".csv";
    ofstream csv(fname);
    for (int j = 0; j < mesh.count; ++j){
        stringstream ss1, ss2;
        auto [distr, width] = get_distribution_per_leaf(N, j);
        ss1 << j + 1 << "; width = " << width << " bottom, ";
        ss2 << j + 1 << "; width = " << width << " top, ";
        for (auto& i: distr){
            double dist = i[1] - i[0];
            if (dist <= 0){
                ss1 << "-, ";
                ss2 << "-, ";
            }
            else {
                ss1 << i[0] << ", ";
                ss2 << i[1] << ", ";
            }
        }
        csv << ss1.str() << "\n" << ss2.str() << endl;
    }
}

void CoaptHeight::compute_leaf_profiles(){
    computeCoaptationField();
    if (m_leaf_profiles_b) return;
    for (int i = 0; i < mesh.count; ++i)
        compute_plane_leaf_coaptation_profile(i);
    m_leaf_profiles_b = true;
}
void CoaptHeight::compute_plane_leaf_coaptation_profile(int mesh_n){
    vector<array<int, 3>> triags;

    auto cn = common_nodes_with_connection(mesh_n);
    auto [pnts1, link1, mapp1] = _get_connected_pnts(m_projection[remap[mesh_n][0]], m_edges[remap[mesh_n][0]]);
    auto [pnts2, link2, mapp2] = _get_connected_pnts(m_projection[remap[mesh_n][1]], m_edges[remap[mesh_n][1]]);
    pnts1 = align_to_dir(pnts1, proj_dir[remap[mesh_n][0]]);
    pnts2 = align_to_dir(pnts2, proj_dir[remap[mesh_n][1]]);
//        auto is_here = [](auto& link, auto& mapp){
//            auto is_here_ = [&link = link, &mapp = mapp](pair<int, int> q){
//                int dd[2] = {mapp[q.first], mapp[q.second]};
//                if (dd[0] > dd[1]) swap(dd[0], dd[1]);
//                return (link.count({dd[0], dd[1]}) > 0);
//            };
//            return is_here_;
//        };
//        string direct = "../Results/res1579026461/";
//        string name;
//        std::ofstream off_file;
//        auto e_to_pnts_ids = [](auto& pnts){
//            auto e_to_pnts1 = [&pnts = pnts](auto& e){
//                array<point_t, 3> arr;
//                for (int i = 0; i < 3; ++i)
//                    arr[i] = VEC(pnts[e[i]][0], pnts[e[i]][1], 0);
//                return arr;
//            };
//            return e_to_pnts1;
//        };

    map<int, int> common;
    for (auto& i: cn){
        common.insert({mapp1[i->contact_elem_id[0]], mapp2[i->contact_elem_id[1]]});
    }
    auto shift = pnts2[(common.begin())->second] - pnts1[(common.begin())->first];
    for(auto& i: pnts2){
        i = i - shift;
    }

    auto p2ar = [](Point_2 p) { return array<double, 2>({p[0], p[1]}); };
    auto n2ar = [&pp = pnts1](pair<int, int> id) {
        Point_2& p = pp[id.first];
        return array<double, 2>({p[0], p[1]});
    };
    auto n2ar2 = [&pp = pnts2](pair<int, int> id) {
        Point_2& p = pp[id.second];
        return array<double, 2>({p[0], p[1]});
    };
//        off_file.open(direct + "test1.off");
//        _savetoOFF(off_file, pnts1, p2ar);
//        off_file.close();
//        off_file.open(direct + "test2.off");
//        _savetoOFF(off_file, pnts2, p2ar);
//        off_file.close();
//        off_file.open(direct + "test3.off");
//        _savetoOFF(off_file, common, n2ar);
//        off_file.close();
//        off_file.open(direct + "test4.off");
//        _savetoOFF(off_file, common, n2ar2);
//        off_file.close();

//        off_file.open(direct + "test1.stl");
//        auto node_conv1 = [](node_t* n) {if (!n->contact_elem_id) return -1; return n->contact_elem_id[0]; };
//        triags = convert_edges_to_elems(node_conv1, is_here(link1, mapp1), mesh_n);
//        _savetoSTL(off_file, triags, e_to_pnts_ids(pnts1));
//        off_file.close();
//
//        auto node_conv2 = [](node_t* n) {if (!n->contact_elem_id) return -1; return n->contact_elem_id[1]; };
//        off_file.open(direct + "test2.stl");
//        triags = convert_edges_to_elems(node_conv2, is_here(link2, mapp2), mesh_n);
//        _savetoSTL(off_file, triags, e_to_pnts_ids(pnts2));
//        off_file.close();


    SewTwoFields sw(pnts1, link1, pnts2, link2, common);
    sw.find_minimum_energy_df(-1, 1.e-4, 1.e-2, 5.e-1, 2000, 2);

//        off_file.open(direct + "test11.stl");
//        auto node_conv11 = [](node_t* n) {if (!n->contact_elem_id) return -1; return n->contact_elem_id[0]; };
//        triags = convert_edges_to_elems(node_conv11, is_here(link1, mapp1), mesh_n);
//        _savetoSTL(off_file, triags, e_to_pnts_ids(pnts1));
//        off_file.close();
//
//        auto node_conv21 = [](node_t* n) {if (!n->contact_elem_id) return -1; return n->contact_elem_id[1]; };
//        off_file.open(direct + "test21.stl");
//        triags = convert_edges_to_elems(node_conv21, is_here(link2, mapp2), mesh_n);
//        _savetoSTL(off_file, triags, e_to_pnts_ids(pnts2));
//        off_file.close();

    auto invmap = [](map<int, int> m){
        map<int, int> res;
        for (const auto& i: m)
            res.insert({i.second, i.first});
        return res;
    };
    auto brd1 = get_minmax_Point_set(pnts1, Vector_2(1, 0));
    auto brd2 = get_minmax_Point_set(pnts2, Vector_2(1, 0));
    if (brd1.first[1] > brd2.first[1]){
        auto inv = [](vector<Point_2>& pp){
            for (auto& i: pp)
                i = Point_2(i[0], -i[1]);
        };
        inv(pnts1);
        inv(pnts2);
    }
    map<int, int> gmap;
    set<pair<int, int>> links;
    auto invmap1 = invmap(mapp1), invmap2 = invmap(mapp2);
    auto invcommon = invmap(common);
    auto& nodes1 = m_nodes[remap[mesh_n][0]], nodes2 = m_nodes[remap[mesh_n][1]];
    set<int> exclude;
    for (const auto& i: common)
        exclude.insert(i.second);
    vector<Point_2> res;
    res.reserve(pnts1.size() + pnts2.size() - cn.size());
    int off = 0;
    for (int i = 0; i < pnts1.size(); ++i) {
        res.push_back(pnts1[i]);
        gmap.insert({off++, nodes1[invmap1[i]]->id});
    }
    for (auto& i: link1){
        array<int, 2> temp = {i.first, i.second};
        if (temp[0] > temp[1]) swap(temp[0], temp[1]);
        links.insert({temp[0], temp[1]});
    }
    map<int, int> loc_ids2;
    for (int i = 0; i < pnts2.size(); ++i){
        if (!exclude.count(i)) {
            res.push_back(pnts2[i]);
            gmap.insert({off++, nodes2[invmap2[i]]->id});
            loc_ids2.insert({i, off-1});
        }
        else loc_ids2.insert({i, invcommon[i]});
    }
    for (auto& i: link2){
        array<int, 2> temp = {loc_ids2[i.first], loc_ids2[i.second]};
        if (temp[0] > temp[1]) swap(temp[0], temp[1]);
        links.insert({temp[0], temp[1]});
    }

    set<pair<int, int>> boundary_links;
    for (const auto& i: m_boundary_edges[remap[mesh_n][0]]){
        array<int, 2> temp = {(int)nodes1[i.first]->id, (int)nodes1[i.second]->id};
        if (temp[0] > temp[1]) swap(temp[0], temp[1]);
        if (links.count({temp[0], temp[1]}))
            boundary_links.insert({temp[0], temp[1]});
    }
    for (const auto& i: m_boundary_edges[remap[mesh_n][1]]){
        array<int, 2> temp = {(int)nodes2[i.first]->id, (int)nodes2[i.second]->id};
        if (temp[0] > temp[1]) swap(temp[0], temp[1]);
        if (links.count({temp[0], temp[1]}))
            boundary_links.insert({temp[0], temp[1]});
    }

//        off_file.open(direct + "test11.off");
//        _savetoOFF(off_file, pnts1, p2ar);
//        off_file.close();
//        off_file.open(direct + "test21.off");
//        _savetoOFF(off_file, pnts2, p2ar);
//        off_file.close();
//        off_file.open(direct + "test31.off");
//        _savetoOFF(off_file, common, n2ar);
//        off_file.close();
//        off_file.open(direct + "test41.off");
//        _savetoOFF(off_file, common, n2ar2);
//        off_file.close();
//        off_file.open(direct + "test_full1.off");
//        _savetoOFF(off_file, res, p2ar);
//        off_file.close();


//        auto check_link = [&links = links, inmp = invmap(gmap)](pair<int, int> q) mutable {
//            int tmp[2] = {inmp[q.first], inmp[q.second]};
//            if (tmp[0] > tmp[1]) swap(tmp[0], tmp[1]);
//            return (links.count({tmp[0], tmp[1]}) > 0);
//        };
//        off_file.open(direct + "test_full1.stl");
//        auto nconv1 = [](node_t* n) {return n->id;};
//        auto inmp = invmap(gmap);
//        auto e_to_pnts3 = [&pnts = res, &inmp = inmp](auto& e){
//            array<point_t, 3> arr;
//            for (int i = 0; i < 3; ++i)
//                arr[i] = VEC(pnts[inmp[e[i]]][0], pnts[inmp[e[i]]][1], 0);
//            return arr;
//        };
//        triags = convert_edges_to_elems(nconv1, check_link, mesh_n);
//        _savetoSTL(off_file, triags, e_to_pnts3);
//        off_file.close();


    m_leaflet_coapt_fields.insert({mesh_n, {move(res), move(links), move(gmap), move(boundary_links)}});
}

set<node_t*> CoaptHeight::common_nodes(int mesh_n){
    int i = mesh_n;
    set<node_t*> common;
    if (mesh.count < 3) {
        common = colData[{i, (i + 1) % 2}];
        return common;
    }
    auto& d1 = colData[remap[i][0]], &d2 = colData[remap[i][1]];
    set_intersection(d1.begin(), d1.end(), d2.begin(), d2.end(), inserter(common, common.begin()));
    return common;
}
set<spring_t*> CoaptHeight::common_edges(int mesh_n){
    set<spring_t*> common;
    auto cn = common_nodes(mesh_n);
    for (int i = 0, cnt = mesh.nets[mesh_n].springs.count; i < cnt; ++i){
        spring_t* spr = mesh.nets[mesh_n].springs.springs[i];
        if (cn.count(spr->ends[0]) && cn.count(spr->ends[1]))
            common.insert(spr);
    }
}
set<node_t*> CoaptHeight::common_nodes_with_connection(int mesh_n){
    auto cn = common_nodes(mesh_n);
    auto &ed1 = m_edges[remap[mesh_n][0]],
            &ed2 = m_edges[remap[mesh_n][1]];
    auto connectivity = [](set<pair<int, int>> s){
        set<int> res;
        for (auto& i: s){
            res.insert(i.first);
            res.insert(i.second);
        }
        return res;
    };
    auto ced1 = connectivity(ed1), ced2 = connectivity(ed2);
    for (auto it = cn.begin(); it != cn.end(); ){
        int i1 = (*it)->contact_elem_id[0], i2 = (*it)->contact_elem_id[1];
        if (!ced1.count(i1) || !ced2.count(i2)){
            it = cn.erase(it);
        }
        else
            ++it;
    }
    return cn;
}
vector<CoaptHeight::Point_2> CoaptHeight::align_to_dir(const vector<Point_2>& pnts, const Vector_2& dir){
    vector<Point_2> aligned;
    aligned.reserve(pnts.size());
    Vector_2 ort(dir[1], -dir[0]);
    for (auto& i: pnts){
        auto v = i - pnts[0];
        aligned.push_back(Point_2(CGAL::scalar_product(ort, v), CGAL::scalar_product(dir, v)));
    }
    return aligned;
}
vector<CoaptHeight::Point_2> CoaptHeight::dealign_from_dir(const vector<Point_2>& aligned, const Vector_2& dir) {
    vector<Point_2> pnts;
    Point_2 z(0, 0);
    pnts.reserve(aligned.size());
    Vector_2 ort(dir[1], -dir[0]);
    for (auto& i: aligned){
        pnts.push_back(z + i[0] * ort + i[1] * dir);
    }
    return pnts;
}

void CoaptHeight::SewTwoFields::set_data(gsl_vector *x){
    array<int, 4> off;
    off[0] = 0;
    off[1] = off[0] + m_p1.size() - m_cn.size();
    off[2] = off[1] + m_cn.size();
    off[3] = off[2] + m_p2.size() - m_cn.size();
    set<int> overlap1, overlap2;
    for (const auto& i: m_cn){
        overlap1.insert(i.first);
        overlap2.insert(i.second);
    }
    int id = 0;
    for (int i = 0; i < m_p1.size(); ++i){
        if (overlap1.count(i)) continue;
        m_map1.insert({i, id / 2});
        gsl_vector_set(x, id++, m_p1[i][0]);
        gsl_vector_set(x, id++, m_p1[i][1]);
    }
    for (auto i: overlap1){
        m_map1.insert({i, id / 2});
        m_map2.insert({m_cn[i], id / 2});
        gsl_vector_set(x, id++, m_p1[i][0]);
        gsl_vector_set(x, id++, m_p1[i][1]);
    }
    for (int i = 0; i < m_p2.size(); ++i){
        if (overlap2.count(i)) continue;
        m_map2.insert({i, id / 2});
        gsl_vector_set(x, id++, m_p2[i][0]);
        gsl_vector_set(x, id++, m_p2[i][1]);
    }
}

void CoaptHeight::SewTwoFields::set_new_points(gsl_vector* x){
    for(int i = 0, cnt = m_p1.size(); i < cnt; ++i) {
        int j = m_map1[i];
        m_p1[i] = Point_2(gsl_vector_get(x, 2 * j), gsl_vector_get(x, 2 * j + 1));
    }
    for(int i = 0, cnt = m_p2.size(); i < cnt; ++i) {
        int j = m_map2[i];
        m_p2[i] = Point_2(gsl_vector_get(x, 2 * j), gsl_vector_get(x, 2 * j + 1));
    }
}

double CoaptHeight::SewTwoFields::compute(const gsl_vector *x){
    double f = 0;
    auto adder = [&f = f, &m_wy = m_wy, &x = x](auto& container, auto& from, auto& _conv) {
        for (const auto &i: container) {
            double cur[4] = {from[i.first][0], from[i.first][1], from[i.second][0], from[i.second][1]};
            int idd[2] = {2*_conv[i.first],2*_conv[i.second]};
            double next[4] = {gsl_vector_get(x, idd[0]), gsl_vector_get(x, idd[0] + 1),
                              gsl_vector_get(x, idd[1]), gsl_vector_get(x, idd[1] + 1)};
            double weights[] = {1, m_wy};
            for (int k = 0; k < 2; ++k) {
                double loc = next[k] - next[k + 2] - (cur[k] - cur[k + 2]);
                f += loc * loc * weights[k];
                //cout << loc * loc * weights[k] << " ";
            }
        }
    };
    adder(m_l1, m_p1, m_map1);
    adder(m_l2, m_p2, m_map2);
    return f;
}

void CoaptHeight::SewTwoFields::compute_df(const gsl_vector *x, gsl_vector *df){
    auto adder = [&](auto& container, auto& from, auto& conv) {
        for (const auto &i: container) {
            int ids[2] = {i.first, i.second};
            int idd[2] = {2*conv[i.first], 2*conv[i.second]};
            double cur[4] = {from[i.first][0], from[i.first][1], from[i.second][0], from[i.second][1]};
            double next[4] = {gsl_vector_get(x, idd[0]), gsl_vector_get(x, idd[0] + 1),
                              gsl_vector_get(x, idd[1]), gsl_vector_get(x, idd[1] + 1)};
            double weights[] = {1, m_wy};
            for (int k = 0; k < 2; ++k) {
                double loc = next[k] - next[k + 2] - (cur[k] - cur[k + 2]);
                double ins = loc * weights[k];
                gsl_vector_set(df, idd[0] + k ,gsl_vector_get(df, idd[0] + k) + ins);
                gsl_vector_set(df, idd[1] + k,gsl_vector_get(df, idd[1] + k) - ins);
            }
        }
    };
    adder(m_l1, m_p1, m_map1);
    adder(m_l2, m_p2, m_map2);
}

double CoaptHeight::SewTwoFields::compute_fdf(const gsl_vector *x, gsl_vector *df){
    double f = 0;
    auto adder = [&](auto& container, auto& from, auto& conv) {
        for (const auto &i: container) {
            int ids[2] = {i.first, i.second};
            int idd[2] = {2*conv[i.first], 2*conv[i.second]};
            double cur[4] = {from[i.first][0], from[i.first][1], from[i.second][0], from[i.second][1]};
            double next[4] = {gsl_vector_get(x, idd[0]), gsl_vector_get(x, idd[0] + 1),
                              gsl_vector_get(x, idd[1]), gsl_vector_get(x, idd[1] + 1)};
            double weights[] = {1, m_wy};
            for (int k = 0; k < 2; ++k) {
                double loc = next[k] - next[k + 2] - (cur[k] - cur[k + 2]);
                double ins = loc * weights[k];
                f += ins * loc;
                gsl_vector_set(df, idd[0] + k ,gsl_vector_get(df, idd[0] + k) + ins);
                gsl_vector_set(df, idd[1] + k,gsl_vector_get(df, idd[1] + k) - ins);
            }
        }
    };
    adder(m_l1, m_p1, m_map1);
    adder(m_l2, m_p2, m_map2);
    return f;
}

void CoaptHeight::SewTwoFields::_fdf4gsl(const gsl_vector *x, void *params, double *f, gsl_vector *df)
{
    gsl_vector_set_all(df, 0);
    SewTwoFields* m = (SewTwoFields*)params;
    *f = m->compute_fdf(x, df);
}
double CoaptHeight::SewTwoFields::_f4gsl (const gsl_vector *x, void *params)
{
    SewTwoFields* m = (SewTwoFields*)params;

    return m->compute(x);
}
void CoaptHeight::SewTwoFields::_df4gsl (const gsl_vector *x, void *params, gsl_vector *df)
{
    gsl_vector_set_all(df, 0);
    SewTwoFields* m = (SewTwoFields*)params;
    m->compute_df(x, df);
}

int CoaptHeight::SewTwoFields::find_minimum_energy_df(int freq, double step_sz, double tol, double epsabs, int maxits, double time){
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

    size_t iter = 0;
    int status;

    const gsl_multimin_fdfminimizer_type *T;
    gsl_multimin_fdfminimizer *s;

    void* par = this;

    int g_N = (m_p1.size() + m_p2.size() - m_cn.size()) * 2;
    gsl_vector *x = gsl_vector_alloc (g_N);
    set_data(x);

    gsl_multimin_function_fdf my_func;

    my_func.n = g_N;
    my_func.f = _f4gsl;
    my_func.df = _df4gsl;
    my_func.fdf = _fdf4gsl;
    my_func.params = par;

    T = gsl_multimin_fdfminimizer_conjugate_fr;
    s = gsl_multimin_fdfminimizer_alloc (T, g_N);

    gsl_multimin_fdfminimizer_set (s, &my_func, x, step_sz, tol);

    _Timer t;
    do
    {
        iter++;
        status = gsl_multimin_fdfminimizer_iterate (s);

        if (status)
            break;

        status = gsl_multimin_test_gradient (s->gradient, epsabs);

        if (status == GSL_SUCCESS && freq >= 0)
            printf ("Minimum found at: %5lu: f() = %10.5f\n", iter, s->f);

        if ((freq > 0) && !(iter % freq))
            printf ("%5lu: f() = %10.5f\n", iter, s->f);

    }
    while (status == GSL_CONTINUE && iter < maxits && t.elapsed() < time);
    if (status != GSL_SUCCESS && freq >= 0)
        printf ("Success didn't reached status = %d: %5lu: f() = %10.5f\n", status, iter, s->f);

    set_new_points(s->x);

    gsl_multimin_fdfminimizer_free (s);
    gsl_vector_free (x);

    return status;
}

tuple<vector<CoaptHeight::Point_2>, set<pair<int, int>>, map<int, int>> CoaptHeight::_get_connected_pnts(const vector<Point_2>& pnts, const set<pair<int, int>>& connect){
    set<int> use_ids;
    for (const auto& i: connect)
        use_ids.insert(i.first), use_ids.insert(i.second);
    vector<Point_2> use_pnts;
    use_pnts.reserve(use_ids.size());
    map<int, int> to_new_ids;
    int counter = 0;
    for (auto i: use_ids) {
        use_pnts.push_back(pnts[i]);
        to_new_ids.insert({i, counter});
        counter++;
    }
    set<pair<int, int>> new_connect;
    for (auto& i: connect){
        new_connect.insert({to_new_ids[i.first], to_new_ids[i.second]});
    }
    return tuple(use_pnts, new_connect, to_new_ids);

}
pair<vector<array<double, 2>>, double> CoaptHeight::get_distribution_per_leaf(int N, int mesh_n){
    compute_leaf_profiles();
    auto& [pnts, links, gmap, borders] = m_leaflet_coapt_fields[mesh_n];
    Vector_2 dir(0, 1), ort(1, 0);
    const double eps = 1e-5;
    double from =  get_minmax_Point_set(pnts, dir).first[1] - eps;
    auto _pp = get_minmax_Point_set(pnts, ort);
    double pp[2] = {_pp.first[0], _pp.second[0]};
    vector<array<double, 2>> distribution;
    distribution.reserve(N+1);
    for (int i = 0; i < N+1; ++i){
        double cf = static_cast<double>(i) / N;
        auto point = from * dir + (pp[0] + cf * (pp[1] - pp[0]))*ort;
        Point_2 l_p(point[0], point[1]);
        auto diam = _clever_field_directed_diam(pnts, links, Line_2(l_p, dir), dir);
        distribution.push_back({diam[0] - from, diam[1] - from});
    }
    return pair(move(distribution), pp[1] - pp[0]);
}
