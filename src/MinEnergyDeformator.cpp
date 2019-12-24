#include "MinEnergyDeformator.h"

class LengthFunctor: public MinEnergyDeformator::EnergyFunctor{
public:
    double weight;
    LengthFunctor (double weigth = 0) : weight{weigth} {}
    virtual double operator()(const net_t& mesh, int edge_id, double df[6]) override {
        double alpha = weight;
        spring_t& sp = *mesh.springs.springs[edge_id];
        double l = LEN(DIF(sp.ends[0]->coord, sp.ends[1]->coord));
        double l_0 = sp.l_0;
        double eps = (l - l_0) / l_0;
        double f = alpha * eps * eps;
        if (df) {
            double dcoef = 2 * alpha * eps / l / l_0;
            for (int i = 0; i < 3; ++i){
                df[i] = (sp.ends[0]->coord.coord[i] - sp.ends[1]->coord.coord[i]) * dcoef;
                df[i+3] = -df[i];
            }
        }
        return f;
    }
};

class DigedralFunctor: public MinEnergyDeformator::EnergyFunctor {
public:
    double weight;
    double scale;
    bool convexity;

    DigedralFunctor (double weight = 0, double scale = 1, bool convexity = true) :
    weight{weight}, scale{scale}, convexity{convexity} {}
    virtual double operator()(const net_t& mesh, int edge_id, double df[12]) override {
        point_t p[4];
        spring_t& sp = *mesh.springs.springs[edge_id];
        p[0] = sp.dihedral[0]->coord;
        p[1] = sp.dihedral[1]->coord;
        p[2] = sp.ends[0]->coord;
        p[3] = sp.ends[1]->coord;
        double* x = df, k;
        if (x) {
            double phi = deriv_f3(p, x, 1.0, scale, &k);
            double weightphik = weight * phi;
            double _edge_coef = -2 * weightphik;
            for (int i = 0; i < 12; ++i)
                x[i] *= _edge_coef;
            return weightphik * phi * k;
        }
        else{
            return weight * get_phi3(p, x, 1.0, scale);
        }
    }

private:
    inline double get_phi(const point_t p[3], point_t res[3], double sc[2]){
        double a2 = SQR_LEN(p[0]);
        double a = sqrt(a2);
        point_t axb = CROSS(p[0], p[1]), cxa = CROSS(p[2], p[0]);
        double abc = DOT(axb, p[2]), axb_axc = -DOT(axb, cxa);
        double phi = abc * a;
        double psi = axb_axc;
        double f = phi / psi;
        sc[0] = phi, sc[1] = -psi;

        return f;
    }
   inline double get_phi2(const point_t p[3], point_t res[3], double beta, double scal){
        double sc[2];
        double x = get_phi(p, res, sc);
        double f = atan2(sc[0], sc[1]) / M_PI;
        double k = ((convexity) ? (f > 0) : (f < 0)) ? 0.1 : (scal * 0.1);

        return k * f * f;
    }

    inline double get_phi3(const point_t p[4], double x[12], double beta, double scal){
        int mask[3] = {3, 0, 1};
        point_t pp[3];
        for (int i = 0; i < 3; ++i)
            pp[i] = ADD(p[mask[i]], -1, p[2]);
        double f = get_phi2(pp, NULL, beta, scal);
        return f;
    }

    double deriv_f1(const point_t p[3], point_t res[3], double sc[2]){
        //f = (a, b, c) * |a| / ([a, b], [a, c])
        //p[0] <-> a
        //p[1] <-> b
        //p[2] <-> c
        double a2 = SQR_LEN(p[0]);
        double a = sqrt(a2);
        double ab = DOT(p[0], p[1]), ac = DOT(p[0], p[2]), bc = DOT(p[1], p[2]);
        point_t axb = CROSS(p[0], p[1]), cxa = CROSS(p[2], p[0]), bxc = CROSS(p[1], p[2]);
        double abc = DOT(axb, p[2]), axb_axc = -DOT(axb, cxa);
        //std::array<point_t&, 3> abc_p = {bxc, cxa, axb};
        double phi = abc * a;
        //point_t abc_adna = SCAL(abc/a, p[0]);
        //for (int i = 1; i < 3; ++i)
        //    SCAL_S(a, &abc_p);
        //point_t phi_da = SCAL_SUM(a, bxc, abc/a, p[0]);
        std::array<point_t, 3> phi_p = {SCAL_SUM(a, bxc, abc/a, p[0]), SCAL(a, cxa), SCAL(a, axb)};
        double psi = axb_axc;
        point_t psi_da;
        for (int i = 0; i < 3; ++i)
            psi_da.coord[i] = 2 * bc * p[0].coord[i] - ac * p[1].coord[i] - ab * p[2].coord[i];
        point_t psi_p[3] = {psi_da, SCAL_SUM(a2, p[2], -ac, p[0]), SCAL_SUM(a2, p[1], -ab, p[0])};
        double f = phi / psi;
        sc[0] = phi, sc[1] = -psi;
        double d1_psi = 1.0/psi, phid2psi = -phi / psi / psi;
        for (int i = 0; i < 3; ++i){
            res[i] = SCAL_SUM(d1_psi, phi_p[i], phid2psi, psi_p[i]);
        }
        return f;
    }

    double deriv_f2(const point_t p[3], point_t res[3], double beta, double scal, double* kk){
        double sc[2];
        double x = deriv_f1(p, res, sc);
        double f = atan2(sc[0], sc[1]) / M_PI;
        double coef_arctg = 1.0 / (1.0 + x * x);
        double k = ((convexity) ? (f > 0) : (f < 0)) ? 0.1 : (scal * 0.1);
        double coef = coef_arctg * k * beta / M_PI;
        for (int i = 0; i < 3; ++i)
            SCAL_S(coef, &res[i]);
        *kk = k;
        return f;
    }

    double deriv_f3(const point_t p[4], double x[12], double beta, double scal, double* kk){
        int mask[3] = {3, 0, 1};
        point_t pp[3];
        for (int i = 0; i < 3; ++i)
            pp[i] = ADD(p[mask[i]], -1, p[2]);
        point_t grad[3];
        double f = deriv_f2(pp, grad, beta, scal, kk);

        if (x != nullptr)
            for (int i = 0; i < 3; ++i) {
                x[3 * 0 + i] = grad[1].coord[i];
                x[3 * 1 + i] = grad[2].coord[i];
                x[3 * 2 + i] = -(grad[0].coord[i] + grad[1].coord[i] + grad[2].coord[i]);
                x[3 * 3 + i] = grad[0].coord[i];;
            }

        return f;
    }

};

class PlaneConstrFunctor: public MinEnergyDeformator::EnergyFunctor {
public:
    point_t normal;
    double m_b;   // (n, x) = b
    double weight;

    PlaneConstrFunctor(double weight, const point_t& normal, double b) :
    weight{weight}, normal{normal}, m_b{b} {}
    virtual double operator()(const net_t& mesh, int node_id, double* df) override {
        double pn = DOT(normal, mesh.vrtx.nodes[node_id]->coord);
        double pnmb = pn - m_b;
        double mpnmb2 = pnmb * 2 * weight;
        bool intersect = (pnmb > 0);
        if (df && intersect){
            for (int i = 0; i < 3; ++i)
                df[i] = mpnmb2 * normal.coord[i];
        }
        return (intersect) ? weight * pnmb * pnmb : 0;
    }
};

void set_default_length_constr(MinEnergyDeformator& m, double weight){
    auto en = std::make_shared<LengthFunctor>(weight);
    m.addEdgeEnergyComponent(en);
//        m.addEdgeEnergyComponent(new LengthFunctor(weight));
}

void set_default_digedral_angle_constr(MinEnergyDeformator& m, double weight, double scale, bool convexity){
    auto en = std::make_shared<DigedralFunctor>(weight, scale, convexity);
    m.addDigedralEnergyComponent(en);
//        m.addDigedralEnergyComponent(new DigedralFunctor(weight, scale, convexity));
}

void set_plane_constr(MinEnergyDeformator& m, double weight, const point_t& normal, double b){
    auto en = std::make_shared<PlaneConstrFunctor>(weight, normal, b);
    m.addNodeEnergyComponent(en);
//        m.addNodeEnergyComponent(new PlaneConstrFunctor(weight, normal, b));
}
