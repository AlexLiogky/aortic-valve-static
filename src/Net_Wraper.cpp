#include "NetCollisionShape.h"
#include "Net_Wraper.h"
#include "computation.h"
#include <array>
#include <vector>
#include <algorithm>
#include <iostream>
#include <assert.h>
#include "NetSpliter.h"

btVector3 p2Bt (const point_t& p){ return btVector3(p.coord[0], p.coord[1],p.coord[2]); }
point_t Bt2p (const btVector3& p){ return point_t_get_point(p.x(), p.y(), p.z()); }

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

NetObject::NetCollisionObject(net_t net, Net_Wraper* body): m_net{net}, m_body{body}
{
    m_nleaf.resize(m_net.vrtx.count);
	for (int i = 0, ni = m_net.vrtx.count; i < ni; ++i)
	{
		node_t* n = m_net.vrtx.nodes[i];
		m_nleaf[i] = m_ndbvt.insert(btDbvtVolume::FromCR(p2Bt(n->coord), m_body->getMargin()), n);
	}

	m_fleaf.resize(m_net.elems.count);
	for (unsigned int i = 0; i < m_net.elems.count; ++i)
	{
		elem_t* f = m_net.elems.elems[i];
		m_fleaf[i] = m_fdbvt.insert(VolumeOf(f, 0), f);
	}

    m_worldTransform.setIdentity();                                                 ////////////neccesary?
	m_collisionShape = new NetCollisionShape(this);
	m_collisionShape->setMargin(m_body->getMargin());
	m_collisionFlags &= ~CF_STATIC_OBJECT;
	updateBounds();
}

void NetObject::updateBounds(){
    if (m_fdbvt.m_root)
	{
		const btVector3& mins = m_fdbvt.m_root->volume.Mins();
		const btVector3& maxs = m_fdbvt.m_root->volume.Maxs();
		const btScalar csm = getCollisionShape()->getMargin();
		const btVector3 mrg = btVector3(csm, csm, csm) * 1;
		m_bounds[0] = mins - mrg;
		m_bounds[1] = maxs + mrg;
	}
	else
	{
		m_bounds[0] = m_bounds[1] = btVector3(0, 0, 0);
	}

}

void NetObject::updateCollisionInfo(){
    ATTRIBUTE_ALIGNED16(btDbvtVolume) vol;

    for (int i = 0, ni = m_net.vrtx.count; i < ni; ++i)
	{
		node_t* n = m_net.vrtx.nodes[i];
		vol = btDbvtVolume::FromCR(p2Bt(n->next), getMargin());
		btDbvtNode* leaf = m_nleaf[i];
		btVector3 velocity = p2Bt(DIF(n->next, n->coord)) * 3;
		btScalar margin = getMargin() * (btScalar)0.25;
		m_ndbvt.update(leaf, vol, velocity, margin);
	}

	#define FN_V(FACE, I) DIF(FACE->vrts[I]->next, FACE->vrts[I]->coord)
	for (int i = 0, fi = m_net.elems.count; i < fi; ++i)
		{
			elem_t* f = m_net.elems.elems[i];
			const point_t v = SCAL(1.0/3, SUM(SUM(FN_V(f, 0), FN_V(f, 1)), FN_V(f, 2)));

			vol = VolumeOf(f, getMargin());
			m_fdbvt.update(m_fleaf[i],
						   vol,
						   p2Bt(v) * 3,
						   getMargin() * (btScalar)0.25);
		}
    #undef FN_V

    m_ndbvt.optimizeIncremental(1);
	m_fdbvt.optimizeIncremental(1);

	updateBounds();
}

void NetObject::CollisionHandler(NetObject* psb){
    if (this->m_body == psb->m_body) return;
    NetObject::Collider docol;
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

void NetObject::CollisionHandler(const btCollisionObject* pco, btVector3* triangle)
{

    NetObject::ColliderStatic docollide;

    ATTRIBUTE_ALIGNED16(btDbvtVolume) volume;
    volume = btDbvtVolume::FromPoints(triangle, 3);
    const btScalar margin = pco->getCollisionShape()->getMargin();
    volume.Expand(btVector3(margin, margin, margin));
    docollide.psb = this;
    for (int i = 0; i < 3; ++i) docollide.m_triangle[i] = triangle[i];
    docollide.mrg = pco->getCollisionShape()->getMargin();

    m_ndbvt.collideTV(m_ndbvt.m_root, volume, docollide);
}

void NetObject::ColliderStatic::Process(const btDbvtNode* leaf)
{
    node_t* node = (node_t*)leaf->data;
    btVector3* bttr = m_triangle;
    point_t n[] = {Bt2p(bttr[0]), Bt2p(bttr[1]), Bt2p(bttr[2])};
    double d = DBL_MAX;
    point_t proj = point_to_triangle_projection(node->coord, n, &d);

    const double m = mrg + LEN(DIF(node->next, node->coord));
    if (d < (m * m))
    {
        //const point_t w = point_t_bary_coord(n[0], n[1], n[2], proj);
        double len = sqrt(d);

        Net_Wraper::RContact c;
        c.m_normal = SCAL(1.0 / len, DIF(node->coord, proj));
        c.m_proj = proj;
        c.m_margin = m;
        c.m_node = node;
        psb->m_body->m_rcontacts.push_back(c);
    }
}

void NetObject::Collider::Process (const btDbvtNode* lnode,
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
            psb[0]->m_body->m_scontacts.push_back(c);
        }
    }
}

void Net_Wraper::splitNet2ColObjs(const std::vector<int>& axisSeq, int depth)
{
    if (depth)
    {
        NetSpliter sp{m_net};
        sp.split(axisSeq, depth);
        std::vector<net_t> preBodies = sp.getBodyPreforms();
        m_objs.reserve(preBodies.size());
        int cnt = 0;
        assert(preBodies[1].vrtx.nodes[0]);
        for (auto& i: preBodies)
            m_objs.push_back(new NetObject(i, this));
    }
    else
    {
        m_objs.push_back(new NetObject(m_net, this));
    }
}

Net_Wraper::Net_Wraper(net_t net, double P, double delta) : m_net{net}, m_P{P}, m_delta{delta}
{
    std::vector<int> axisSeq({1, 1});
    splitNet2ColObjs(axisSeq, 0);
}

void Net_Wraper::addObjs2World(btCollisionWorld* wrld)
{
    for (int i = 0, bc = m_objs.size(); i < bc; ++i)
        wrld->addCollisionObject(m_objs[i], btBroadphaseProxy::DefaultFilter, btBroadphaseProxy::AllFilter);
}

void Net_Wraper::updateCollisionInfo()
{
    for (int i = 0, bc = m_objs.size(); i < bc; ++i)
        m_objs[i]->updateCollisionInfo();

    m_scontacts.resize(0);
    m_rcontacts.resize(0);
}

double Net_Wraper::computeFreeNexts(){
    if (net_is_static(this->m_net)) return 0;
    return compute_free_nexts(m_net, m_P, m_delta);
}

point_t BaryEval(point_t a, point_t b , point_t c, point_t w)
{
    point_t res;
    for (int i = 0; i < 3; ++i)
        res.coord[i] = a.coord[i]*w.coord[0] + b.coord[i]*w.coord[1] + c.coord[i]*w.coord[2];
    return res;
}

void Net_Wraper::solveRContacts(Net_Wraper* psb)
{
    for (int i = 0, ni = psb->m_rcontacts.size(); i < ni; ++i)
        {
            const RContact& c = psb->m_rcontacts[i];
            const point_t& nr = c.m_normal;
            const point_t& p = c.m_proj;
            node_t& n = *c.m_node;

            const point_t vr = DIF(n.next, n.coord);
            point_t corr = get_zero_point();
            double dot = DOT(vr, nr);
            if (dot < 0)
            {
                const double j = c.m_margin - (DOT(nr, n.next) - DOT(nr, p));
                ADD_S(&corr, j, c.m_normal);
            }
            //corr -= ProjectOnPlane(vr, nr) * c.m_friction;
            ADD_S(&n.next, 1, corr);
        }
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
                const double j = c.m_margin - (DOT(nr, n.next) - DOT(nr, p));
                ADD_S(&corr, j, c.m_normal);
            }
            //corr -= ProjectOnPlane(vr, nr) * c.m_friction;
            ADD_S(&n.next, c.m_cfm[0], corr);
            ADD_S(&f.vrts[0]->next, -c.m_cfm[1] * c.m_weights.coord[0], corr);
            ADD_S(&f.vrts[1]->next, -c.m_cfm[1] * c.m_weights.coord[1], corr);
            ADD_S(&f.vrts[2]->next, -c.m_cfm[1] * c.m_weights.coord[2], corr);
        }
}




