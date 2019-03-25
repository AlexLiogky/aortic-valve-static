#ifndef NET_WRAPER_H
#define NET_WRAPER_H

#include "nets.h"
#include "BulletCollision/BroadphaseCollision/btDbvt.h"
#include "LinearMath/btAlignedObjectArray.h"

class Net_Wraper
{

public:
    struct Collider : public btDbvt::ICollide
    {
        Net_Wraper* psb[2];
        double mrg;

        virtual void Process (const btDbvtNode* lnode,
                const btDbvtNode* lface) override;
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

typedef btAlignedObjectArray<btDbvtNode*> tDbvtArray;
typedef btAlignedObjectArray<SContact> tSContactArray;

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
        double getMargin() { return m_mrg; }
        double computeFreeNexts();
        void updateCollisionInfo();
        static void solveSContacts(Net_Wraper* body);


private:
    tSContactArray m_scontacts;        // Soft contacts
    tDbvtArray m_nleaf;
    tDbvtArray m_fleaf;
    btDbvt m_ndbvt;                    // Nodes tree
	btDbvt m_fdbvt;                    // Faces tree
	net_t m_net;
	double m_mrg = 0.2;
	double m_P;
	double m_delta;
};

btVector3 p2Bt (point_t p);

#endif // NET_WRAPER_H
