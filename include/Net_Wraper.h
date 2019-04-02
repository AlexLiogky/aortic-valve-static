#ifndef NET_WRAPER_H
#define NET_WRAPER_H

#include "nets.h"
#include "BulletCollision/BroadphaseCollision/btDbvt.h"
#include "LinearMath/btAlignedObjectArray.h"
#include "BulletCollision/CollisionDispatch/btCollisionObject.h"
#include "BulletCollision/CollisionDispatch/btCollisionWorld.h"
#include <vector>
#include "BulletSoftBody/btSparseSDF.h"

btVector3 p2Bt (const point_t& p);
point_t Bt2p (const btVector3& p);

class Net_Wraper
{

public:
    class NetCollisionObject: public btCollisionObject
    {
    public:
        struct Collider : public btDbvt::ICollide
        {
            NetCollisionObject* psb[2];
            double mrg;

            virtual void Process (const btDbvtNode* lnode,
                    const btDbvtNode* lface) override;
        };


        struct ColliderStatic : public btDbvt::ICollide
        {
            NetCollisionObject* psb;
            btVector3 m_triangle[3];
            double mrg;

            virtual void Process (const btDbvtNode* leaf) override;
        };
        typedef btAlignedObjectArray<btDbvtNode*> tDbvtArray;


    public:
        Net_Wraper* m_body;
        net_t m_net;
        tDbvtArray m_nleaf;
        tDbvtArray m_fleaf;
        btDbvt m_ndbvt;                    // Nodes tree
        btDbvt m_fdbvt;                    // Faces tree
        btVector3 m_bounds[2];

        NetCollisionObject(net_t net, Net_Wraper* backPtr);
        void CollisionHandler(NetCollisionObject* psb);
        void CollisionHandler(const btCollisionObject* pco, btVector3* triangle);
        double getMargin() { return m_body->getMargin(); }
        void updateCollisionInfo();
        virtual void getAabb(btVector3& aabbMin, btVector3& aabbMax) const
        {
            aabbMin = m_bounds[0];
            aabbMax = m_bounds[1];
        }
        void updateBounds();
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

typedef btAlignedObjectArray<SContact> tSContactArray;
typedef btAlignedObjectArray<RContact> tRContactArray;
typedef btAlignedObjectArray<NetCollisionObject*> ObjArray;

public:
        Net_Wraper(net_t net, double P, double delta);

        double getMargin() { return m_mrg; }
        double computeFreeNexts();
        void addObjs2World(btCollisionWorld* wrld);
        void updateCollisionInfo();
        static void solveSContacts(Net_Wraper* body);
        static void solveRContacts(Net_Wraper* body);
private:
        void splitNet2ColObjs(const std::vector<int>& axisSeq, int depth);

private:
    ObjArray m_objs;
    tSContactArray m_scontacts;        // Soft contacts
    tRContactArray m_rcontacts;
	net_t m_net;
	double m_mrg = 0.2;
	double m_P;
	double m_delta;
};

using NetObject = Net_Wraper::NetCollisionObject;

#endif // NET_WRAPER_H
