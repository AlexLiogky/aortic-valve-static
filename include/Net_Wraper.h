#ifndef NET_WRAPER_H
#define NET_WRAPER_H

#include "nets.h"
#include "BulletCollision/BroadphaseCollision/btDbvt.h"
#include "LinearMath/btAlignedObjectArray.h"
#include "BulletCollision/CollisionDispatch/btCollisionObject.h"
#include "BulletSoftBody/btSparseSDF.h"

btVector3 p2Bt (const point_t& p);
point_t Bt2p (const btVector3& p);

class Net_Wraper : public btCollisionObject
{

public:
    struct Collider : public btDbvt::ICollide
    {
        Net_Wraper* psb[2];
        double mrg;

        virtual void Process (const btDbvtNode* lnode,
                const btDbvtNode* lface) override;
    };


    struct ColliderStatic : public btDbvt::ICollide
    {
        Net_Wraper* psb;
        btVector3 m_triangle[3];
        double mrg;

        virtual void Process (const btDbvtNode* leaf) override;
    };

struct SContact
	{
		node_t* m_node;         // Node
		elem_t* m_face;         // Face
		point_t m_weights;      // Weigths
		point_t m_normal;       // Normal
		double m_margin;        // Margin
		double m_cfm[2];        // Constraint force mixing
	};

struct RContact
	{
		node_t* m_node;         // Node
		point_t m_proj;         // Projection
		point_t m_normal;       // Normal
		double m_margin;        // Margin
	};

typedef btAlignedObjectArray<btDbvtNode*> tDbvtArray;
typedef btAlignedObjectArray<SContact> tSContactArray;
typedef btAlignedObjectArray<RContact> tRContactArray;

public:
        Net_Wraper(net_t net, double P, double delta);
        /*Net_Wraper(Net_Wraper&& net){                                    /////////////////////////change after/////////
            m_scontacts = net.m_scontacts;
            m_nleaf = std::move(net.m_nleaf);
            m_fleaf = std::move(net.m_fleaf);
            m_ndbvt = std::move(net.m_ndbvt);
            m_fdbvt = std::move(net.m_fdbvt);
            m_net = net.m_net;
            m_mrg = net.m_mrg;
            m_P = net.m_P;
            m_delta = net.m_delta;
        }*/
        void CollisionHandler(Net_Wraper* psb);
        void CollisionHandler(const btCollisionObject* pco, btVector3* triangle);
        double getMargin() { return m_mrg; }
        double computeFreeNexts();
        void updateCollisionInfo();
        static void solveSContacts(Net_Wraper* body);
        static void solveRContacts(Net_Wraper* body);
        virtual void getAabb(btVector3& aabbMin, btVector3& aabbMax) const
        {
            aabbMin = m_bounds[0];
            aabbMax = m_bounds[1];
        }
        void updateBounds();

public:
    btSparseSdf<3> m_sparsesdf;
private:
    tSContactArray m_scontacts;        // Soft contacts
    tRContactArray m_rcontacts;
    tDbvtArray m_nleaf;
    tDbvtArray m_fleaf;
    btDbvt m_ndbvt;                    // Nodes tree
	btDbvt m_fdbvt;                    // Faces tree
	net_t m_net;
	double m_mrg = 0.2;
	double m_P;
	double m_delta;
public:
	btVector3 m_bounds[2];
};

#endif // NET_WRAPER_H
