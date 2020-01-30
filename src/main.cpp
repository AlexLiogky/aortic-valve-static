#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <math.h>
#include <string.h>
#include <sys/time.h>

#include <SDL/SDL.h>
#include <GL/gl.h>
#include <GL/glu.h>

#include "nets.h"
#include "stress.h"
#include "precomputation.h"
#include "computation.h"
#include "postcomputation.h"
#include "intersection.h"
#include "save-data.h"
#include "format_in.h"
#include "separate.h"
#include "sew_leaf.h"
#include "bound-box.h"
#include "camera.h"
//#include "Solver.h"
#include "World.h"
#include "InputProcessor.h"
#include "Meassure.h"
#include <string>
#include <sstream>
#include <fstream>
#include <iostream>
#include <map>
#include <algorithm>

#include "MinEnergyDeformator.h"
using namespace std;

struct MyLog{
    void open( const std::string &filename,
               ios_base::openmode mode = ios_base::out ){
        if (activated)
            log.open(filename, mode);
    }
    void close(){
        if (activated)
            log.close();
    }
    void set_activation(bool state) { activated = state; }
    template <typename T>
    MyLog& operator<< (const T& value){
        if (activated)
            log << value;
        return *this;
    }
private:
    std::ofstream log;
    bool activated = true;
};

InputProcessor gPrms;
MyLog gLog;

//######################################################################
long get_msec_time(struct timeval start, struct timeval end){
	long seconds  = end.tv_sec  - start.tv_sec;
    long useconds = end.tv_usec - start.tv_usec;
	long mtime = ((seconds) * 1000 + useconds / 1000.0) + 0.5;
	return mtime;
}

long double get_msec_time_(struct timeval start, struct timeval end){
	long double seconds  = end.tv_sec  - start.tv_sec;
    long double useconds = end.tv_usec - start.tv_usec;
	long double mtime = (seconds) * 1000 + useconds / 1000.0;
	return mtime;
}

void print_max_min(nets_t nets){
	int nets_c = nets.count, i;
	double max = 0, min = nets.nets[0].springs.springs[0][0].l_0;
	for(i = 0; i < nets_c; i++){
		net_t net = nets.nets[i];
		int s_c = net.springs.count, j;
		for (j = 0; j < s_c; j++){
			double len = net.springs.springs[j][0].l_0;
			if (max < len) max = len;
			if (min > len) min = len;
		}
	}
	printf("max-spr = %lg, min-spr = %lg\n", max, min);
}

void save_data(nets_t nets, char* file_name){
	char name[1024] = {0};
	strcpy(name, file_name);
	int len = strlen(file_name);
	nets_t net = nets_t_get_net(1);
	for (unsigned char i = 0; i < nets.count; i++){
		name[len] = '0' + i;
		net.nets[0] = nets.nets[(int)i];
		printf("name = <%s>\n", name);
		to_stl(net, name);
	}
}

nets_t download_aorta(const char* file_name){
	nets_t aorta = download_nets_from_file(file_name);
	precomputation(aorta);
	return aorta;
}

point_t _triag_center(point_t triag[3]){
    point_t a[3];
    for (int i = 0; i < 3; ++i)
        a[i] = DIF(triag[(i+1) %3], triag[(i+2) %3]);

    double s2 = SQR_LEN(CROSS(a[0], a[1])); //4s^2
    point_t res = ZERO();
    for (int i = 0; i < 3; ++i)
        ADD_S(&res, -SQR_LEN(a[i]) / s2 / 2 * DOT(a[(i + 1) %3], a[(i + 2) % 3]), triag[i]);
    return res;
}

string to_string(point_t p){
    return "(" + to_string(p.coord[0]) + ", " + to_string(p.coord[1]) + ", " + to_string(p.coord[2]) + ")";
}

void min_energ_to_bnd(net_t& leaf, point_t att[2], point_t blood, point_t center){
    point_t n[2];
    for (int i = 0; i < 2; ++i)
        n[i] = NORM(CROSS(blood, DIF(center, att[i])));
    SCAL_S(-1, &n[1]);

    MinEnergyDeformator m(leaf);
    SewEnergyParams& sep = gPrms.sep;
    set_plane_constr(m, sep.plane_w, n[0], DOT(center, n[0]));
        gLog << "  plane: n = " << to_string(n[0]) << ", x0 = " << to_string(center) << "\n";
    set_plane_constr(m, sep.plane_w, n[1], DOT(center, n[1]));
        gLog << "  plane: n = " << to_string(n[1]) << ", x0 = " << to_string(center) << "\n";
    set_default_length_constr(m, sep.sqr_length_w);
    set_default_digedral_angle_constr(m, /*sep.digedral_angle_w*/0.7, /*sep.convexity_w*/1.0);
    point_t wind = NORM(SCAL_SUM(-0.3, NORM(blood), 1.0, SCAL_SUM(0.5, n[0], 0.5, n[1])));
        gLog << "  plane: wind = " << to_string(wind) << "\n";
    set_isotrop_force(m, gPrms.sep.force_w, wind);
    SewEnergyParams::EnergyMinimizerParams& emp = sep.sp;
    m.find_minimum_energy_df(emp.freq, emp.step_sz, emp.tol, emp.epsabs, emp.maxits, emp.time);
}

void shift_free_mesh(net_t leaf, double scale, const point_t shift){
    for (int i = 0; i < leaf.vrtx.count; ++i)
    {
        if (is_fix(leaf.vrtx.nodes[i]->state))
            continue;
        ADD_S(&leaf.vrtx.nodes[i]->coord, scale, shift);
    }
}

int fit_mid_points(bnds_t& bnds, point_t cntr, point_t axis, point_t mids[3][2], double scale){
    point_t e0 = NORM(CROSS(axis, CROSS(DIF(bnds.bnds[0].line[0], cntr), axis)));
    point_t e1 = NORM(CROSS(axis, e0));
    double phi[3][2] = {};
    double dphi = 0;
    for (int ii = 0 ; ii < bnds.cnt; ++ii) {
        point_t e00 = NORM(CROSS(axis, CROSS(DIF(bnds.bnds[ii].line[0], cntr), axis)));
        point_t e10 = NORM(CROSS(axis, e00));
        point_t loc = DIF(bnds.bnds[ii].line[0], cntr);
        phi[ii][0] = phi[ii][1] = std::atan2(DOT(e10, loc), DOT(e00, loc));
        for (int i = 1; i < bnds.bnds[ii].pnt_cnt; ++i) {
            loc = DIF(bnds.bnds[ii].line[i], cntr);
            double _phi = std::atan2(DOT(e10, loc), DOT(e00, loc));
            if (_phi < phi[ii][0]) phi[ii][0] = _phi;
            else if (_phi > phi[ii][1]) phi[ii][1] = _phi;
        }
        if (phi[ii][1] - phi[ii][0] > M_PI) {
            phi[ii][0] += 2 * M_PI;
            std::swap(phi[ii][0], phi[ii][1]);
        }
        double ddphi = atan2(DOT(e00, e1), DOT(e00, e0));
        phi[ii][0] += ddphi, phi[ii][1] += ddphi;
        if (phi[ii][0] < dphi) dphi = phi[ii][0];
    }
    for (int i = 0; i < 6; ++i)
        phi[i%3][i%2] -= dphi;
    std::vector<std::pair<double*, int>> s = {{phi[0], 0}, {phi[1], 1}, {phi[2], 2}};
    struct {
        bool operator()(std::pair<double*, int> a, std::pair<double*, int> b) const
        {
            return a.first[0] < b.first[0];
        }
    } customLess;
    std::sort(s.begin(), s.end(), customLess);
    int nmap[3] = {0, 2, 1};    //сумма индексов кусков минус один
    double pmids[3];
    double pdifs[3];
    int err = 0;

        pmids[nmap[s[0].second + s[1].second - 1]] = (s[0].first[1] + s[1].first[0]) / 2;
        pdifs[nmap[s[0].second + s[1].second - 1]] = fabs(s[0].first[1] - s[1].first[0]) / 2;

    if (s[0].first[1] > s[1].first[0]) {
        std::cout << "Warning: fields " << s[0].second << ", " << s[1].second << " doesn't separable\n", err = -1;
        pdifs[nmap[s[0].second + s[1].second - 1]] = 0;
    }
    {
        pmids[nmap[s[1].second + s[2].second - 1]] = (s[1].first[1] + s[2].first[0]) / 2;
        pdifs[nmap[s[1].second + s[2].second - 1]] = fabs(s[1].first[1] - s[2].first[0]) / 2;
    }
    if (s[1].first[1] > s[2].first[0]) {
        std::cout << "Warning: fields " << s[1].second << ", " << s[2].second << " doesn't separable\n", err = -1;
        pdifs[nmap[s[1].second + s[2].second - 1]] = 0;
    }
    {
        pmids[nmap[s[2].second + s[0].second - 1]] = (s[0].first[0] -( 2*M_PI - s[2].first[1])) / 2;
        pdifs[nmap[s[2].second + s[0].second - 1]] = (s[0].first[0] +( 2*M_PI - s[2].first[1])) / 2;
    }
    if (s[2].first[1] > 2 * M_PI){
        std::cout << "Warning: fields " << s[2].second << ", " << s[0].second << " doesn't separable\n", err = -1;
        pdifs[nmap[s[2].second + s[0].second - 1]] = 0;
    }
    for (int i = 0; i < 3; ++i)
        pmids[i] += dphi;
    point_t n;
    for (int i = 0; i < 6; ++i) {
        double psi = pmids[(i % 3 + i % 2 + 2) % 3] -
                     (1 - 2 * (i % 2)) * pdifs[(i % 3 + i % 2 + 2) % 3] / 2.0;
        mids[i % 3][(i + 1) % 2] = SUM(cntr,  SCAL_SUM(scale * cos(psi), e0, scale * sin(psi), e1));
//        n = CROSS(axis, SCAL_SUM(scale * cos(psi), e0, scale * sin(psi), e1));
//        std::cout << to_string(n) << std::endl;
    }

    return err;
}

void  energetical_sew_leaf(net_t leafs[3], point_t att[3][2], point_t blood, bnds_t* bnds = nullptr){
    int neigh_id[6] = {0, 1, 2, 3, 4, 5};
    double lens[6];
    for (int i = 0; i < 6; ++i)
    for (int j = 0; j < 6; ++j){
        if (i == j) continue;
        double len2 = SQR_LEN(DIF(att[i/2][i%2], att[j/2][j%2]));
        if (neigh_id[i] == i || lens[i] > len2){
            lens[i] = len2;
            neigh_id[i] = j;
        }
    } //for (int i = 0; i < 6; ++i) lens[i] = sqrt(lens[i]);
    point_t mids[6];
    for (int i = 0; i < 6; ++i)
        mids[i] = SCAL_SUM(0.5, att[i/2][i%2], 0.5, att[neigh_id[i]/2][neigh_id[i]%2]);

    point_t triag[3];
    triag[0] = mids[0];
    int notes[6] = {0, neigh_id[0], -1, -1, -1, -1};
        for (int i = 1, j = 1; i < 6 && j < 3; ++i){
            int flag = 1;
            for (int k = 0; k < 2 * j; ++k)
                if (notes[k] == i) {
                    flag = 0;
                    break;
                }
            if (!flag) continue;
            triag[j++] = mids[i];
            notes[2 * j] = i;
            notes[2 * j + 1] = neigh_id[i];
        }

    point_t center = _triag_center(triag);
    double R = LEN(DIF(triag[0], center));

    for (int i = 0; i < 3; ++i)
        shift_free_mesh(leafs[i], /*-1.0*/gPrms.sep.shift * R, NORM(blood));

    /*nets_t leaflet = nets_t_get_net(3);
    for (int i = 0; i < 3; ++i)
        leaflet.nets[i] = leafs[i];
    leaflet.count = 3;
	to_stl(leaflet, "shift");
	exit(-2);*/

    gLog << "Sewing started: \n";
    if (bnds != nullptr)
        fit_mid_points(*bnds, center, blood, att, R);

    for (int i = 0; i < 3; ++i){
            gLog << " Leaflet["<< i << "]: ";
        point_t shift = ZERO();//SUM(DIF(att[i][0], mids[i * 2]), DIF(att[i][1], mids[i * 2 + 1]));
            gLog << "shift = " << to_string(shift) << "\n";
        min_energ_to_bnd(leafs[i], att[i], blood, SUM(center, shift));
    }

}

nets_t get_system(nets_t aorta, const char* bnd, point_t blood_direction, const char* leaf1, const char* leaf2, const char* leaf3){
	point_t shift = blood_direction;

	point_t init_points[2], final_points[3][2];
	point_t init_shift = point_t_get_point(0, -9, 0);
	point_t fiber_direction = point_t_get_point(0, 1, 0);

	point_t commissur[3];
	bnds_t bnds = read_bnds(bnd);
	for (int i = 0; i < bnds.cnt; i++)
		commissur[i] = bnds.bnds[i].line[0];
	/*point_t bottom[3] = {{{-11.941, -36.7436, 2.81578}}, {{-3.33692, -17.0099, 14.5751}}, {{-14.1942, -21.4188, -4.19314}}};
	sew_line_t* sews = init_sew_line_on_aorta(aorta, commissur, bottom, 3);
	//sew_line_t_plate_dump(&sews[0]);
	point_t p = {{-7.4382, -28.391, 20.944}};
	sew_line_t_add_point(p, &sews[1]);
	sew_line_t_add_point(VEC(-4.994, -26.491, 18.874), &sews[1]);
	sew_line_t_add_point(VEC(-3.515, -23.861, 18.608), &sews[1]);
	sew_line_t_add_point(VEC(-16.265, -13.279, 10.909), &sews[1]);
	bnds_t newbnds = {3, calloc(3, sizeof(line_t))};
	for (int i = 0; i < 3; i++)
		newbnds.bnds[i] = sew_line_t_get_line_on_aorta(&sews[i]);
	printf("count points in leaf[0] = %d\n", newbnds.bnds[1].pnt_cnt);
	//newbnds.bnds[0] = sew_line_t_get_line_on_aorta(&sews[0]);
	line_t_dump(newbnds.bnds[1]);
	bnds = newbnds;*/


	nets_t leaflet1 = formated_in(leaf1);
	nets_t_set_relax_state(leaflet1, fiber_direction);
	to_stl(leaflet1, (gPrms.ro.res_dir + "relaxed.stl").c_str());
	init_points[0] = leaflet1.nets[0].vrtx.nodes[1][0].coord;
	init_points[1] = leaflet1.nets[0].vrtx.nodes[4][0].coord;
	final_points[0][0] = commissur[1];
	final_points[0][1] = commissur[0];
	sew_leaflet_to_aorta(leaflet1.nets[0], init_points, init_shift,\
						aorta.nets[0], final_points[0], shift, 4, bnds.bnds[0]);

	nets_t leaflet2 = formated_in(leaf2);
	nets_t_set_relax_state(leaflet2, fiber_direction);
	init_points[0] = leaflet2.nets[0].vrtx.nodes[1][0].coord;
	init_points[1] = leaflet2.nets[0].vrtx.nodes[4][0].coord;
	final_points[1][0] = commissur[2];
	final_points[1][1] = commissur[1];
	sew_leaflet_to_aorta(leaflet2.nets[0], init_points, init_shift,\
						aorta.nets[0], final_points[1], shift, 4, bnds.bnds[1]);

	nets_t leaflet3 = formated_in(leaf3);
	nets_t_set_relax_state(leaflet3, fiber_direction);
	init_points[0] = leaflet3.nets[0].vrtx.nodes[1][0].coord;
	init_points[1] = leaflet3.nets[0].vrtx.nodes[4][0].coord;
	final_points[2][0] = commissur[0];
	final_points[2][1] = commissur[2];
	sew_leaflet_to_aorta(leaflet3.nets[0], init_points, init_shift,\
						aorta.nets[0], final_points[2], shift, 4, bnds.bnds[2]);

	int n_objects = 4;
	nets_t leaflet = nets_t_get_net(n_objects);
	leaflet.nets[0] = leaflet1.nets[0];
	leaflet.nets[1] = leaflet2.nets[0];
	leaflet.nets[2] = leaflet3.nets[0];
	leaflet.nets[3] = aorta.nets[0];
	for (int i = 0; i < gPrms.ots.size() && 3; ++i)
	    net_t_set_thickness(leaflet.nets[i], gPrms.ots[i].thickness);
	net_t_set_state(&leaflet.nets[3], 1);
	precomputation(leaflet);

	leaflet.count = 3;
	to_stl(leaflet, (gPrms.ro.res_dir + "start").c_str());
	leaflet.count = n_objects;
	energetical_sew_leaf(leaflet.nets, final_points, blood_direction, &bnds);
	leaflet.count = 3;
	to_stl(leaflet, (gPrms.ro.res_dir + "minim").c_str());
	leaflet.count = n_objects;
	//exit(-1);

	return leaflet;
	//return aorta;
}

nets_t get_valve_from_system(nets_t system){
	nets_t valve = nets_t_get_net(system.count - 1);
	for (unsigned int i = 0; i < system.count - 1; i++)
		valve.nets[i] = system.nets[i];

	return valve;
}

void print_nets_statistic(nets_t leaflet, point_t shift){
    const double recognition = 1.05 * 4;
                                        set_contact_recognition_resolution(recognition * 0.1);
	printf("Contact_Resolution = %lg\n", get_Contact_Resolution());
        gLog << "Contact Resolution = " << get_Contact_Resolution()<< "\n";
	nets_t curleaflet = leaflet;//get_valve_from_system(leaflet);
	print_statistic(curleaflet);
	double l_free[10];
	for (unsigned int i = 0; i < leaflet.count; ++i)
        l_free[i] = net_t_get_len_free_edge(leaflet.nets[i]);

	double h = nets_t_get_coapt_depth(curleaflet, &shift);
	int _pow = -1;
	                                    set_contact_recognition_resolution(recognition * 0.35);
	double h_c = nets_t_get_coapt_intersect_depth(curleaflet, &_pow);
                                        set_contact_recognition_resolution(recognition * 0.1);
    printf("l_free[leaf] = {%lg, %lg, %lg}\n", l_free[0], l_free[1], l_free[2]);
	printf("%s: h = %lg, h_c = %lg\n", gPrms.ro.res_dir.c_str(), h, h_c);
        gLog << "l_free[leaf] = { " <<  l_free[0] << ", " << l_free[1] << ", " << l_free[2] << " }\n";
        gLog << "0-th level:\n" << " h = " << h << "\n";
        gLog << " h_c = " << h_c << ", count of central points = " << _pow << "\n";

	for (int i = 0; i < 1; i++){
		curleaflet = create_next_hierarchical_nets(curleaflet);
            gLog << i + 1 << "-th level:\n";
		//auto_set_contact_recognition_consts(curleaflet, NULL);

		                                    set_contact_recognition_resolution(recognition * 0.35);
		double h_c = nets_t_get_coapt_intersect_depth(curleaflet, &_pow);
		                                    set_contact_recognition_resolution(recognition * 0.1);
		double h = nets_t_get_coapt_depth(curleaflet, &shift);
            gLog << " h = " << h << "\n";
            gLog << " h_c = " << h_c << ", count of central points = " << _pow << "\n";
		double f_area[3] = {};
		double S_max = 0, S_min = 0;
		printf("Compute full leaflet coaptation area:\n");
            gLog << " Full leaflet coaptation area:\n";
		for (int i = 0; i < 3; i++){
			 f_area[i] = net_t_get_coapt_area(curleaflet.nets[i]);
			 if (i == 0) S_max = f_area[i], S_min = f_area[i];
			 S_max = (S_max < f_area[i]) ? f_area[i] : S_max;
			 S_min = (S_min > f_area[i]) ? f_area[i] : S_min;
			 printf("  full area[%d] = %lg\n", i, f_area[i]);
                gLog << "  full area[" << i << "] = " << f_area[i] << "\n";
		 }
		double p_area[6] = {};
		printf("Compute partition leaflet coaptation area:\n");
            gLog << " Partition leaflet coaptation area:\n";
		for (int i = 0; i < 3; i++){
			p_area[2 * i] = nets_t_get_coapt_area_from_to(curleaflet, i%3, (i+1)%3);
			p_area[2 * i + 1] = nets_t_get_coapt_area_from_to(curleaflet, (i+1)%3, i%3);
			printf("partion area[%d -> %d] = %lg\n", i%3, (i+1)%3, p_area[2 * i]);
			printf("partion area[%d -> %d] = %lg\n", (i + 1)%3, (i)%3, p_area[2 * i + 1]);
                gLog << "  partion area[" << i%3 << " -> " << (i+1)%3 << "] = " << p_area[2 * i] << "\n";
                gLog << "  partion area[" << (i + 1)%3 << " -> " << (i)%3 << "] = " << p_area[2 * i + 1] << "\n";
		}
        double h_mid[10];
        for (unsigned int i = 0; i < leaflet.count; ++i)
            h_mid[i] = (p_area[2 * ((i + 2) % 3) + 0] + p_area[2 * i + 1]) / l_free[i];
                gLog << " h_mid[leaf] = { " <<  h_mid[0] << ", " << h_mid[1] << ", " << h_mid[2] << " }\n";
		printf("%s: h = %lg, h_c = %lg, S_max = %lg, S_min = %lg\n", gPrms.ro.res_dir.c_str(), h, h_c, S_max, S_min);
		printf("%s: h_mid[leaf] = { %lg, %lg, %lg }\n", gPrms.ro.res_dir.c_str(), h_mid[0], h_mid[1], h_mid[2]);
	}
}

point_t double_to_rgb(double x, double scale)
{
    x /= scale;
    point_t p = {};
    if (x < 0)
    {
        p.coord[2] = 1;
        return p;
    }
    if (x < 0.5)
    {
        p.coord[2] = 1 - 2* x;
        p.coord[1] = 2 * x;
        return p;
    }
    if (x < 1)
    {
        p.coord[1] = 1 - 2* (x - 0.5);
        p.coord[0] = 2 * (x - 0.5);
        return p;
    }
    p.coord[0] = 1;
    return p;
}

point_t get_color(const net_t& net, node_t* node)
{
    double r = 0;
    for (int i = 0, cnt = node->cnt_springs; i < cnt; ++i)
    {
        spring_t& spr = *net.springs.springs[node->springs_id[i]];
        r += fabs(spr.l - spr.l_0) / spr.l_0;
    }
    r /= node->cnt_springs;
    const double max_deform = 0.6;
    return double_to_rgb(r, max_deform);
}

void renderSoftBody(const net_t net, point_t fcolor){
    glColor3d(fcolor.coord[0], fcolor.coord[1], fcolor.coord[2]);
    glBegin(GL_TRIANGLES);
    for (unsigned int i = 0; i < net.elems.count; ++i)
    {
        elem_t* elem = net.elems.elems[i];
        for(int j = 0; j < 3; ++j)
        {
            fcolor = get_color(net, elem->vrts[j]);
            glColor3d(fcolor.coord[0], fcolor.coord[1], fcolor.coord[2]);
            point_t x = elem->vrts[j]->coord;
            glVertex3d(x.coord[0], x.coord[1], x.coord[2]);
        }
    }
    glEnd();

    glColor3f(0.6, 0.6, 0.6);
    glBegin(GL_LINES);
    for (unsigned int i = 0; i < net.springs.count; ++i)
    {
        spring_t* spr = net.springs.springs[i];
        for(int j = 0; j < 2; ++j)
        {
            point_t x = spr->ends[j]->coord;
            glVertex3d(x.coord[0], x.coord[1], x.coord[2]);
        }
    }
    glEnd();
}


void display(world_t* world, camera_t* cam)
{
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
	glLoadIdentity();
	camera_t_Control(cam);
	//drawSkybox(50);
	camera_t_UpdateCamera(cam);
	for (unsigned i = 0; i < world->dynamic_nets.count; ++i)
        renderSoftBody(world->dynamic_nets.nets[i], point_t_get_point(0.296, 0.221, 0.231));

    for (unsigned i = 0; i < world->static_nets.count; ++i)
        renderSoftBody(world->static_nets.nets[i], point_t_get_point(205.0/256, 200.0/256, 239.0/256));

}

template <class TT>
void displayT(TT* s, camera_t* cam)
{
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
	glLoadIdentity();
	camera_t_Control(cam);
	//drawSkybox(50);
	camera_t_UpdateCamera(cam);
	nets_t& nets = s->getDynamicNets();
	for (unsigned i = 0; i < nets.count; ++i)
        renderSoftBody(nets.nets[i], point_t_get_point(0.296, 0.221, 0.231));

    nets_t& nets1 = s->getStaticNets();
    for (unsigned i = 0; i < nets1.count; ++i)
        renderSoftBody(nets1.nets[i], point_t_get_point(205.0/256, 200.0/256, 239.0/256));

}

void run_comupation(world_t* world){
    SDL_Init(SDL_INIT_EVERYTHING);
	SDL_SetVideoMode(1280,800,32,SDL_OPENGL);
	Uint32 _start;
	SDL_Event event;
	int running=1;
	float angle=45;
	glClearColor(0,0,0,1);
	glMatrixMode(GL_PROJECTION);
		glLoadIdentity();
		gluPerspective(angle,1280.0/800.0,1,1000);
	glMatrixMode(GL_MODELVIEW);
	glEnable(GL_DEPTH_TEST);
	camera_t* cam = camera_t_construct(point_t_get_point(0.800262, -2.22123, 3.87111));//-45.8, -22, 34.6));
	cam->camPitch = 33.2;//-8;
	cam->camYaw = 9;//308;

	while(running)
	{
		_start=SDL_GetTicks();
		while(SDL_PollEvent(&event))
		{
			switch(event.type)
			{
				case SDL_QUIT:
					running=0;
					break;
				case SDL_KEYDOWN:
					switch(event.key.keysym.sym)
					{
						case SDLK_ESCAPE:
							running=0;
							break;
						case SDLK_y:
							camera_t_mouseIn(cam, 0);
							break;
						case SDLK_SPACE:
                        {
                            printf("current possition: "); point_t_dump(cam->loc); printf("\n");
                            printf("camPitch = %lg, camYaw = %lg\n", cam->camPitch, cam->camYaw);
							break;
                        }
                        default: break;
					}
					break;
				case SDL_MOUSEBUTTONDOWN:
				{
                    camera_t_mouseIn(cam, !camera_t_isMouseIn(cam));
					break;
                }
                default: break;

			}
		}
		compute_nets_time(1000.0/60, world, 1000);
        display(world, cam);
		SDL_GL_SwapBuffers();
		if(1000.0/60>SDL_GetTicks()-_start)
			SDL_Delay(1000.0/60-(SDL_GetTicks()-_start));
	}

	SDL_Quit();
	camera_t_destruct(cam);
}

template<class TT>
void run_comupationT(TT* s, unsigned int maxtime_sec = 3600 * 24){
    SDL_Init(SDL_INIT_EVERYTHING);
	SDL_SetVideoMode(1280,800,32,SDL_OPENGL);
	Uint32 _start, start;
	SDL_Event event;
	int running=1;
	int compute=1;
	//int render_mode=1;
	float angle=45;
	glClearColor(0,0,0,1);
	glMatrixMode(GL_PROJECTION);
		glLoadIdentity();
		gluPerspective(angle,1280.0/800.0,1,1000);
	glMatrixMode(GL_MODELVIEW);
	glEnable(GL_DEPTH_TEST);
	auto& orig = gPrms.cp.origin;
	camera_t* cam = camera_t_construct(point_t_get_point(orig[0], orig[1], orig[2]));//-45.8, -22, 34.6));//3.35825, 9.48692, -91.156));//(-45.8, -22, 34.6));//2.79514, 3.93931, 33.3655));//-3.68549, 11.6567, -2.46957));//-45.8, -22, 34.6));
	cam->camPitch = gPrms.cp.pitch;//1;//-8;//4.4;//-62.8;//-8;
	cam->camYaw = gPrms.cp.yaw;//318.6;//308;//7.4;//178;//308;
	//point_t_dump(camera_t_getVector(cam));

	start = SDL_GetTicks();
	unsigned int maxtime_msec = maxtime_sec * 1000;
	while(running && SDL_GetTicks() - start < maxtime_msec)
	{
		_start=SDL_GetTicks();
		while(SDL_PollEvent(&event))
		{
			switch(event.type)
			{
				case SDL_QUIT:
					running=0;
					break;
				case SDL_KEYDOWN:
					switch(event.key.keysym.sym)
					{
						case SDLK_ESCAPE:
							running=0;
							break;
						case SDLK_y:
							camera_t_mouseIn(cam, 0);
							break;
						case SDLK_SPACE:
                        {
                            printf("current possition: "); point_t_dump(cam->loc); printf("\n");
                            printf("camPitch = %lg, camYaw = %lg\n", cam->camPitch, cam->camYaw);
                            compute ^= 1;
							break;
                        }
                        case SDLK_EQUALS:
                        {
                            solver_t data = s->getSolverData();
                            data.delta /= 0.9;
                            s->setSolverData(data);
                            printf("delta = %lg\n", data.delta);
                                gLog << "new delta = " << data.delta << "\n";
                            break;
                        }
                        case SDLK_MINUS:
                        {
                            solver_t data = s->getSolverData();
                            data.delta *= 0.9;
                            s->setSolverData(data);
                            printf("delta = %lg\n", data.delta);
                                gLog << "new delta = " << data.delta << "\n";
                            break;
                        }
                        case SDLK_s:
                        {
                            to_stl(s->m_dynamic_nets, (gPrms.ro.res_dir + "cur_net").c_str());
                            break;
                        }
                        case SDLK_BACKSPACE:
                            exit(-1);
                        default: break;
					}
					break;
				case SDL_MOUSEBUTTONDOWN:
				{
                    camera_t_mouseIn(cam, !camera_t_isMouseIn(cam));
					break;
                }
                default: break;

			}
		}
		double dt = 1000.0/60;
		if (compute)
            s->compute_nets_time(dt, 1000);
        displayT<TT>(s, cam);
		SDL_GL_SwapBuffers();
		if(dt>SDL_GetTicks()-_start)
			SDL_Delay(dt-(SDL_GetTicks()-_start));
	}

	SDL_Quit();
	camera_t_destruct(cam);
}

void prepare_anisotrop_three_leaflet(nets_t valve, point_t direction)
{
    int reserve = 0;
    for (int i = 0, n_cnt = valve.count; i < n_cnt; ++i)
        reserve += valve.nets[i].vrtx.count;
    std::vector<point_t> p_storage; p_storage.reserve(reserve);
    for (int i = 0, n_cnt = valve.count; i < n_cnt; ++i)
    for (int j = 0, v_cnt = valve.nets[i].vrtx.count; j < v_cnt; ++j)
    {
        point_t p = valve.nets[i].vrtx.nodes[j]->coord;
        p_storage.push_back(p);
        double R = sqrt(p.coord[0] * p.coord[0] + p.coord[1] * p.coord[1]);
//        double sin_phi = p.coord[0] / R, cos_phi = p.coord[1] / R;
//        assert(fabs(sin_phi) <= 1 && fabs(cos_phi) <= 1);
        double phi = atan2(p.coord[1], p.coord[0]);
        p.coord[0] = R * phi;
        p.coord[1] = 0;
        valve.nets[i].vrtx.nodes[j]->coord = p;
    }

    precomputation(valve);
    nets_t_set_relax_state(valve, direction);

    for (int i = 0, k = 0, n_cnt = valve.count; i < n_cnt; ++i)
    for (int j = 0, v_cnt = valve.nets[i].vrtx.count; j < v_cnt; ++j)
        valve.nets[i].vrtx.nodes[j]->coord = p_storage[k++];
}

//######################################################################
//######################################################################
int main(int argc, char* argv[]){
//    std::cout << "Start" << std::endl;
//    nets_t res = download_nets_from_file("/home/alex/Desktop/MV/valve model-static/Aortic-Valve/Results/res1579025760/result.nts");
//    to_stl(res, "/home/alex/Desktop/MV/valve model-static/Aortic-Valve/Results/res1579025760/downloaded");
//    return 0;

    gPrms.InputProcessorInit(argc, argv);
    return testMeassure(gPrms);
    if (!gPrms.ro.use) gLog.set_activation(false);
        gLog.open(gPrms.ro.log_name, ios::trunc);
        gLog << to_string(gPrms);


    nets_t aorta = download_aorta(gPrms.a_in.aorta_file.c_str());
	printf("elems.count = %u\n", aorta.nets[0].elems.count);
	printf("Contact_Resolution = %lg\n", get_Contact_Resolution());
	auto& bl_flow = gPrms.a_in.main_flow_direction;
	point_t shift = point_t_get_point(bl_flow[0], bl_flow[1], bl_flow[2]);//(0.714, 0.127, -0.688);//(5, 1, -3.3);//(7.77, 1.73, -6.05);//(5, 1, -3.3);
	nets_t leaflet = get_system(aorta, gPrms.a_in.leaflet_boundaries.c_str(), shift,
	        gPrms.l_ins[0].leaflet_file.c_str(), gPrms.l_ins[1].leaflet_file.c_str(), gPrms.l_ins[2].leaflet_file.c_str());
	nets_t test = nets_t_get_net(leaflet.count-1);
    for (int i = 0; i < leaflet.count-1; ++i) test.nets[i] = leaflet.nets[i];
    if (gPrms.ro.use)
        to_stl(test, (gPrms.ro.res_dir + "minim").c_str());

        gLog << "Elems per net = " << test.nets[0].elems.count << "\n";
        gLog << "Nodes per net = " << test.nets[0].vrtx.count << "\n";
        gLog << "Springs per net = " << test.nets[0].springs.count << "\n";
        gLog << "Initial full area per net = " << net_t_get_full_area(test.nets[0]) << "\n";
    nets_t static_nets = aorta;
    nets_t dynamic_nets = test;

    printf("len = %lg\n", net_t_get_len_free_edge(test.nets[0]));
        gLog << "Initial free edge len per net = " << net_t_get_len_free_edge(test.nets[0]) << "\n";
        gLog << "Fixed len net = " << net_t_get_len_fix_edge(test.nets[0]) << "\n";
//TODO: сейчас модель хранится как параметр задачи, но нужно сделать так, чтобы модель была параметром объекта и притом имела память о своих параметрах
    solver_t solver_data = {gPrms.sp.delta, gPrms.sp.eps, gPrms.ots[0].elastic_model_type}; //1e-7//25e-7
        gLog << "eps = " << solver_data.eps << "\n";
        gLog << "ElasticType = " << solver_data.ElasticModelType << "\n";
        gLog << "delta = " << solver_data.delta << "\n";
    wrld_cnd_t conditions = {gPrms.ots[0].pressure};
        gLog << "P = " << conditions.P << "\n";
    double coef1 = gPrms.sp.max_possible_recommend_shift_scale, coef2 = gPrms.sp.max_possible_shift_scale;
//TODO: DONE: сделать приём массива марджинов для всех тел в World
    World s(dynamic_nets, static_nets, conditions, solver_data, coef2*Allow_shift, coef1*Max_shift, gPrms.ots[0].collision_margin);
    struct timeval start1, end1;
    gettimeofday(&start1, NULL);
    run_comupationT(&s, gPrms.sp.max_time);
    gettimeofday(&end1, NULL);
    printf("Time of computation = %ld ms\n", get_msec_time(start1, end1));
    printf("Time of elastic computation = %lg ms\n", gt_elastic);
        gLog << "Final full area per net = " << net_t_get_full_area(test.nets[0]) << "\n";
        gLog << "Final free edge len per net = " << net_t_get_len_free_edge(test.nets[0]) << "\n";
        gLog << "Time of computation = " << get_msec_time(start1, end1) << "\n";
        gLog << "Time of elastic computation = " << gt_elastic << "\n";

    gPrms.ro.test_data = gPrms.ro.res_dir + "result";   //TODO: to configfile
    if (gPrms.ro.use)
        save_nets_to_file(s.getDynamicNets(), gPrms.ro.test_data.c_str());

//    World::ColissionType coaptor = s.getCollision(0.11);
//    for (auto& i: coaptor){
//        std::cout << std::get<0>(i) << " " << std::get<1>(i) << " ";
//        point_t_dump(std::get<2>(i)->coord);
//    }
//    return 0;
    /*std::map<int, node_t*> coapt;
    for (auto& i: coaptor)
        if (!coapt.count(i->id)) coapt.insert(std::pair<int, node_t*>(i->id, i));
    std::ofstream off_file(gPrms.ro.res_dir + "coapt1"+ "-" + gPrms.ro.postfix +".off");
    off_file << "OFF\n\n";
    off_file << coapt.size() << " 0 0\n";
    for (auto& j: coapt)
    {
        node_t* i = j.second;
        off_file << i->coord.coord[0] << " " << i->coord.coord[1] << " " << i->coord.coord[2] << "\n";
    }
    off_file.close();
    std::ofstream connect(gPrms.ro.res_dir + "map1" + "-" + gPrms.ro.postfix);
    for (auto& j: coapt)
    {
        node_t* i = j.second;
        connect << i->id << " " << i->coord.coord[0] << " " << i->coord.coord[1] << " " << i->coord.coord[2] << "\n";
    }*/
//TODO: сделать сохранялку которая определяет тип сохраняемого файла по названию
    if (gPrms.ro.use) {
        if (gPrms.ro.divide_leaflets)
            for (int i = 0; i < std::min((int) gPrms.ro.leaf_names.size(), (int) dynamic_nets.count); ++i)
                to_stl1(dynamic_nets.nets[i], gPrms.ro.leaf_names[i].c_str());
        else
            to_stl(dynamic_nets, gPrms.ro.leaf_names[0].c_str());
    }

    print_nets_statistic(dynamic_nets, /*point_t_get_point(0, 0, 1)*/shift);

    return 0;

    nets_t dynamic_nets1 = create_next_hierarchical_nets(dynamic_nets);
    precomputation(dynamic_nets1);
    dynamic_nets1 = create_next_hierarchical_nets(dynamic_nets1);
    precomputation(dynamic_nets1);
    World s1(dynamic_nets1, static_nets, conditions, solver_data, coef2*Allow_shift, coef1*Max_shift, gPrms.ots[0].collision_margin);
    gettimeofday(&start1, NULL);
    run_comupationT(&s, gPrms.sp.max_time);
    gettimeofday(&end1, NULL);
    printf("Time of computation = %ld ms\n", get_msec_time(start1, end1));
    printf("Time of elastic computation = %lg ms\n", gt_elastic);
//    gLog << "Final full area per net = " << net_t_get_full_area(test.nets[0]) << "\n";
//    gLog << "Final free edge len per net = " << net_t_get_len_free_edge(test.nets[0]) << endl;
    gLog << "Time of computation = " << get_msec_time(start1, end1) << "\n";
    gLog << "Time of elastic computation = " << gt_elastic << "\n";

//    if (gPrms.ro.divide_leaflets)
//        for (int i = 0; i < std::min((int)gPrms.ro.leaf_names.size(), (int)dynamic_nets.count); ++i)
//            to_stl1(dynamic_nets1.nets[i], gPrms.ro.leaf_names[i].c_str());
//    else
//        to_stl(dynamic_nets1, gPrms.ro.leaf_names[0].c_str());
    print_nets_statistic(dynamic_nets1, /*point_t_get_point(0, 0, 1)*/shift);


    return 0;
}
//######################################################################
//######################################################################
/*TODO:
 * т.к. теперь есть exact_box, то можно создать замену для функции box_t_local_node_to_net_projection(...)
 * сделать проверку правильности последовательности (по state) в compute_nets(...)
 */

/*  TODO:
 * добавить поле-аккумулятор приложенной силы к узлу
 * вынести вычисление силы упругости как силы вычисляемой по данным всей сети
 * вынести полезные предвычисления индуцируемые той или иной силой упругости в отдельный массив по сетке
 */

 /* можно отключить update_node у неподвижных вершин */
 /*центральная кооптация считается неверно, нужно понять, почему*/

 /*изменена контактная сила*/

/*
 * сделать полный рефакторинг кода
 * */

/*
 * имеет смысл делать константы сил контакта
 * пропорциональными sqrt(S) ~ contact_resolution
 * k = k * get_contact_resolution()
 */

/*nets_t leaflet = formated_in(INPUT);
	nets_t_set_relax_state(leaflet, point_t_get_point(0, 0, 1));
	point_t shift = point_t_get_point(0, 0, 1);

	double P = 1; //mm Hg
	double delta = 3e-7;//3e-7;
	printf("delta = %lg\n", delta);
	int max_nsteps = 3000, freq = 150;
	double eps = 0;
	//compute_nets(leaflet, P, delta, max_nsteps, eps, freq);
	precomputation(leaflet);
	P = 80, max_nsteps = 24000;
	compute_nets(leaflet, P, delta, max_nsteps, eps, freq);
	leaflet = create_next_hierarchical_nets(leaflet);
	leaflet = create_next_hierarchical_nets(leaflet);
	leaflet = create_next_hierarchical_nets(leaflet);

	printf("%s: coapt_depth = %lg\n", RES_STL, nets_t_get_coapt_depth(leaflet, &shift));

	to_stl(leaflet, "result");*/


