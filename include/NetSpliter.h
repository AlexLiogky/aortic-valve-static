#ifndef NETSPLITER_H_INCLUDED
#define NETSPLITER_H_INCLUDED

#include <array>
#include <vector>
#include <map>
#include <stdlib.h>
#include <iostream>

class NetSpliter
{
public:
    struct SpFace
    {
        std::array<int, 3> face;
        std::array<bool, 3> own;
        elem_t* backPtr;
    };
    struct SpNode
    {
        point_t node;
        node_t* backPtr;
    };
    struct SpBounds
    {
        int f_st; //face first index
        int f_nd; //index of last element + 1
        int n_st;
        int n_nd;
    };

    using FaceArr = std::vector<SpFace*>;
    using NodeArr = std::vector<SpNode*>;

private:
    std::vector<int> m_FaceArrOffs;
    std::vector<bool> m_f_in;
    FaceArr m_faces;
    std::vector<int> m_NodeArrOffs;
    NodeArr m_nodes;

    template <unsigned char Axis>
    static bool comp(SpNode* n1, SpNode* n2)
        { return n1->node.coord[Axis] < n2->node.coord[Axis]; }

    //using CmpNodeFnc = bool(*)(SpNode*, SpNode*);

    static bool (*get_comp(int axis)) (SpNode*, SpNode*)
    {
        switch(axis)
        {
            case 0: return comp<0>;
            case 1: return comp<1>;
            case 2: return comp<2>;
            default: return comp<0>;
        }
    }

public:
    NetSpliter(net_t& net)
    {
        node_t** ns = net.vrtx.nodes;
        elem_t** es = net.elems.elems;

        m_nodes.reserve(net.vrtx.count);
        for (int i = 0, nc = net.vrtx.count; i < nc; ++i)
            m_nodes.push_back(new SpNode{ns[i]->coord, ns[i]});

        m_faces.reserve(net.elems.count);
        m_f_in.reserve(net.elems.count);
        for (int i = 0, fc = net.elems.count; i < fc; ++i)
        {
            m_faces.push_back(new SpFace{std::array<int, 3>({(int)es[i]->vrts[0]->id, (int)es[i]->vrts[1]->id, (int)es[i]->vrts[2]->id}), std::array<bool, 3>({true, true, true}), es[i]});
            m_f_in.push_back(true);
            /*std::array<3, point_t*>* face = new std::array<3, point_t*>;
            for (int j = 0; j < 3; ++j)
                face->[j] = &(es[i]->vrts[j]->coord);
            faces.push_back(face);*/
        }
    }

    void split(const std::vector<int>& axisSeq, int depth)
    {
        SpBounds b = {0, (int)m_faces.size(), 0, (int)m_nodes.size()};
        m_NodeArrOffs.push_back(0);
        m_NodeArrOffs.push_back(m_nodes.size());
        m_FaceArrOffs.push_back(0);
        m_FaceArrOffs.push_back(m_faces.size());

        _split(b, axisSeq, depth, 0);
    }

    std::vector<net_t> getBodyPreforms(){
        std::sort(m_NodeArrOffs.begin(), m_NodeArrOffs.end());
        std::sort(m_FaceArrOffs.begin(), m_FaceArrOffs.end());
        std::vector<net_t> bodies;
        bodies.reserve(m_NodeArrOffs.size());
        for (int i = 0, nc = m_NodeArrOffs.size() - 1; i < nc; ++i)
        {
            elems_t e = elems_t_construct(m_FaceArrOffs[i+1] - m_FaceArrOffs[i]);
            for (int j = m_FaceArrOffs[i]; j < m_FaceArrOffs[i+1]; ++j)
                e.elems[j - m_FaceArrOffs[i]] = m_faces[j]->backPtr;
            vrtx_t n = vrtx_t_construct(m_NodeArrOffs[i+1] - m_NodeArrOffs[i]);
            for (int j = m_NodeArrOffs[i]; j < m_NodeArrOffs[i+1]; ++j)
                n.nodes[j - m_NodeArrOffs[i]] = m_nodes[j]->backPtr;

            springs_t doomy = springs_t_construct(0);
            //net_t net = net_t_get(n, e, doomy);
            bodies.push_back(net_t_get(n, e, doomy));
        }


        return std::move(bodies);
    }

    ~NetSpliter()
    {
        for (auto i: m_nodes)
            delete i;
        for (auto i: m_faces)
            delete i;
    }

private:
    void _split(SpBounds bnds, const std::vector<int>& axisSeq, int depth, int step)
    {
        if (depth <= 0 || step < 0 || !axisSeq.size() || ! ((bnds.n_nd - bnds.n_st) / 2) || !((bnds.f_nd - bnds.f_st) / 2)) return;

        int axis = axisSeq[step % axisSeq.size()];
        std::sort(&m_nodes[bnds.n_st], &m_nodes[bnds.n_nd], get_comp(axis));
        int mid = (bnds.n_st + bnds.n_nd) / 2;
        /*std::set<int> n_half;
        for (int i = bnds.n_st; i < mid; ++i)
            n_half.insert(m_nodes[i]->backPtr->id);*/

        for (int i = bnds.f_st; i < bnds.f_nd; ++i)
        {
            int own = 0, in_half = 0;
            for (int j = 0; j < 3; ++j)
            {
                if (m_faces[i]->own[j])
                {
                    ++own;
                    in_half += (m_faces[i]->face[j] < mid);
                }
            }
            if (own)
            {
                if (own - in_half <= own / 3)
                {
                    for (int j = 0; j < 3; ++j)
                        m_faces[i]->own[j] &= (m_faces[i]->face[j] < mid);
                    m_f_in[i] = true;
                }
                else
                {
                    for (int j = 0; j < 3; ++j)
                        m_faces[i]->own[j] &= !(m_faces[i]->face[j] < mid);
                    m_f_in[i] = false;
                }
            }
        }
        int jjj = bnds.f_st;
        for (int i = bnds.f_st; i < bnds.f_nd; ++i)
        {
            if (m_f_in[i] && i != jjj) std::swap<SpFace*>(m_faces[i], m_faces[jjj]);
            if (m_f_in[i]) ++jjj;
        }

        m_NodeArrOffs.push_back(mid);
        m_FaceArrOffs.push_back(jjj);

        SpBounds b1 = {bnds.f_st, jjj, bnds.n_st, mid};
        SpBounds b2 = {jjj, bnds.f_nd, mid, bnds.n_nd};
        _split(b1, axisSeq, depth - (step == (int)axisSeq.size() - 1), (step + 1) % axisSeq.size());
        _split(b2, axisSeq, depth - (step == (int)axisSeq.size() - 1), (step + 1) % axisSeq.size());
    }
public:
    static net_t getExtendedNodeArray(const net_t net)
    {
        std::map<node_t*, int> nd_id;
        for (int i = 0, j = 0; i < (int)net.elems.count * 3; ++i)
        {
            node_t* node = net.elems.elems[i/3]->vrts[i%3];
            j += nd_id.insert(std::make_pair(node, j)).second;
        }

        vrtx_t n = vrtx_t_construct(nd_id.size());
        for (auto& i: nd_id)
        {
            const node_t node = *i.first;
            const double* p = node.coord.coord;
            n.nodes[i.second] = node_t_construct(p[0], p[1], p[2], 0.5, node.state, i.second);
            n.nodes[i.second]->id = i.second;
        }

        elems_t e = elems_t_construct(net.elems.count);
        for (int i = 0; i < (int)net.elems.count; ++i)
        {
            int id0 = nd_id.find(net.elems.elems[i]->vrts[0])->second;
            int id1 = nd_id.find(net.elems.elems[i]->vrts[1])->second;
            int id2 = nd_id.find(net.elems.elems[i]->vrts[2])->second;
            e.elems[i] = elem_t_construct(n.nodes[id0], n.nodes[id1], n.nodes[id2], i);
            e.elems[i]->coef = net.elems.elems[i]->coef;
            e.elems[i]->cntr_mass = net.elems.elems[i]->cntr_mass;
        }
        free(net.vrtx.nodes);

        return net_t_get(n, e, springs_t_construct(0));
    }

};

#endif // NETSPLITER_H_INCLUDED
