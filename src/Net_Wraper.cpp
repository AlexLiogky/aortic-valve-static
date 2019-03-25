#include "Net_Wraper.h"
#include "computation.h"

btVector3 p2Bt (point_t p){ return btVector3(p.coord[0], p.coord[1],p.coord[2]); }

static inline btDbvtVolume VolumeOf(const elem_t* f,
									double margin)
{
	const btVector3 tmp[] = {p2Bt(f->vrts[0]->coord),
                             p2Bt(f->vrts[1]->coord),
                             p2Bt(f->vrts[2]->coord)};
	const btVector3* pts[] = {&tmp[0], &tmp[1], &tmp[2]};
	btDbvtVolume vol = btDbvtVolume::FromPoints(pts, 3);
	vol.Expand(btVector3(margin, margin, margin));
	return (vol);
}

Net_Wraper::Net_Wraper(net_t net, double P, double delta) : m_net{net}, m_P{P}, m_delta{delta}
{

	m_nleaf.resize(m_net.vrtx.count);
	for (int i = 0, ni = m_net.vrtx.count; i < ni; ++i)
	{
		node_t* n = m_net.vrtx.nodes[i];
		m_nleaf[i] = m_ndbvt.insert(btDbvtVolume::FromCR(p2Bt(n->coord), m_mrg), n);
	}

	m_fleaf.resize(m_net.elems.count);
	for (int i = 0; i < m_net.elems.count; ++i)
	{
		elem_t* f = m_net.elems.elems[i];
		m_fleaf[i] = m_fdbvt.insert(VolumeOf(f, 0), f);
	}
}

double Net_Wraper::computeFreeNexts(){
    if (net_is_static(this->m_net)) return 0;
    return compute_free_nexts(m_net, m_P, m_delta);
}

void Net_Wraper::updateCollisionInfo(){
    ATTRIBUTE_ALIGNED16(btDbvtVolume) vol;

    for (int i = 0, ni = m_net.vrtx.count; i < ni; ++i)
	{
		node_t* n = m_net.vrtx.nodes[i];
		vol = btDbvtVolume::FromCR(p2Bt(n->next), m_mrg);
		btDbvtNode* leaf = m_nleaf[i];
		btVector3 velocity = p2Bt(DIF(n->next, n->coord)) * 3;
		btScalar margin = m_mrg * (btScalar)0.25;
		m_ndbvt.update(leaf, vol, velocity, margin);
	}

	#define FN_V(FACE, I) DIF(FACE->vrts[I]->next, FACE->vrts[I]->coord)
	for (int i = 0, fi = m_net.elems.count; i < fi; ++i)
		{
			elem_t* f = m_net.elems.elems[i];
			const point_t v = SCAL(1.0/3, SUM(SUM(FN_V(f, 0), FN_V(f, 1)), FN_V(f, 2)));

			vol = VolumeOf(f, m_mrg);
			m_fdbvt.update(m_fleaf[i],
						   vol,
						   p2Bt(v) * 3,
						   m_mrg * (btScalar)0.25);
		}
    #undef FN_V

    m_scontacts.resize(0);

    m_ndbvt.optimizeIncremental(1);
	m_fdbvt.optimizeIncremental(1);

}

void Net_Wraper::CollisionHandler(Net_Wraper* psb){
    if (this == psb) return;
    Net_Wraper::Collider docol;
    docol.mrg = getMargin() + psb->getMargin();
    /* psb0 nodes vs psb1 faces	*/
    if (!net_is_static(this->m_net)){
    docol.psb[0] = this;
    docol.psb[1] = psb;
    docol.psb[0]->m_ndbvt.collideTT(docol.psb[0]->m_ndbvt.m_root,
                                    docol.psb[1]->m_fdbvt.m_root,
                                    docol);
    }
    /* psb1 nodes vs psb0 faces	*/
    if (!net_is_static(psb->m_net)){
    docol.psb[0] = psb;
    docol.psb[1] = this;
    docol.psb[0]->m_ndbvt.collideTT(docol.psb[0]->m_ndbvt.m_root,
                                    docol.psb[1]->m_fdbvt.m_root,
                                    docol);
    }
}

void Net_Wraper::Collider::Process (const btDbvtNode* lnode,
                                    const btDbvtNode* lface)
{
    node_t* node = (node_t*)lnode->data;
    elem_t* face = (elem_t*)lface->data;
    double d = DBL_MAX;
    point_t proj = node_to_elem_projection(node, face, &d);

    const double m = mrg + LEN(DIF(node->next, node->coord)) * 2;
    if (d < (m * m))
    {
        const point_t n[] = {face->vrts[0]->coord, face->vrts[1]->coord, face->vrts[2]->coord};
        const point_t w = point_t_bary_coord(n[0], n[1], n[2], proj);
        const double ma = 1;
        btScalar mb = 1;
        if (net_is_static(psb[1]->m_net))
        {
            mb = 0;
        }
        const double ms = ma + mb;
        if (ms > 0)
        {
            Net_Wraper::SContact c;
            c.m_normal = SCAL(1.0 / sqrt(d), DIF(node->coord, proj));
            c.m_margin = m;
            c.m_node = node;
            c.m_face = face;
            c.m_weights = w;
            c.m_cfm[0] = ma / ms * 1;
            c.m_cfm[1] = mb / ms * 1;
            psb[0]->m_scontacts.push_back(c);
        }
    }
}

point_t BaryEval(point_t a, point_t b , point_t c, point_t w)
{
    point_t res;
    for (int i = 0; i < 3; ++i)
        res.coord[i] = a.coord[i]*w.coord[0] + b.coord[i]*w.coord[1] + c.coord[i]*w.coord[2];
    return res;
}

void Net_Wraper::solveSContacts(Net_Wraper* psb)
{
    for (int i = 0, ni = psb->m_scontacts.size(); i < ni; ++i)
        {
            const SContact& c = psb->m_scontacts[i];
            const point_t& nr = c.m_normal;
            node_t& n = *c.m_node;
            elem_t& f = *c.m_face;
            const point_t p = BaryEval(f.vrts[0]->next,
                                         f.vrts[1]->next,
                                         f.vrts[2]->next,
                                         c.m_weights);
            const point_t q = BaryEval(f.vrts[0]->coord,
                                        f.vrts[1]->coord,
                                         f.vrts[2]->coord,
                                         c.m_weights);
            const point_t vr = DIF(DIF(n.next, n.coord), DIF(p, q));
            point_t corr = get_zero_point();
            double dot = DOT(vr, nr);
            if (dot < 0)
            {
                const double j = c.m_margin - (DOT(nr, n.coord) - DOT(nr, p));
                ADD_S(&corr, j, c.m_normal);
            }
            //corr -= ProjectOnPlane(vr, nr) * c.m_friction;
            ADD_S(&n.next, c.m_cfm[0], corr);
            ADD_S(&f.vrts[0]->next, -c.m_cfm[1] * c.m_weights.coord[0], corr);
            ADD_S(&f.vrts[1]->next, -c.m_cfm[1] * c.m_weights.coord[1], corr);
            ADD_S(&f.vrts[2]->next, -c.m_cfm[1] * c.m_weights.coord[2], corr);
        }
}


