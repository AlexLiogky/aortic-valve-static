#ifndef NET_WRAPER_H
#define NET_WRAPER_H

#include "nets.h"
#include "BulletCollision/BroadphaseCollision/btDbvt.h"
#include "LinearMath/btAlignedObjectArray.h"
#include "BulletCollision/CollisionDispatch/btCollisionObject.h"
#include "BulletCollision/CollisionDispatch/btCollisionWorld.h"
#include <vector>
#include "BulletSoftBody/btSparseSDF.h"
#include <utility>
#include <tuple>

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
    private:
        std::vector<std::pair<node_t*, Net_Wraper*>> _m_data1;
        std::vector<std::pair<elem_t*, Net_Wraper*>> _m_data2;
    };

    struct SContact
	{
		node_t* m_node;         // Node
		elem_t* m_face;         // Face
		point_t m_weights;      // Weigths
		point_t m_normal;       // Normal
		double m_margin;        // Margin
		double m_cfm[2];        // Constraint force mixing
        Net_Wraper* m_nobstr;   // Node body
        Net_Wraper* m_fobstr;   // Face body
	};

    struct RContact
	{
		node_t* m_node;         // Node
		point_t m_proj;         // Projection
		point_t m_normal;       // Normal
		point_t m_tnormal;      // Normal of triangle
		double m_margin;        // Margin
        Net_Wraper* m_obstr;    // Obstruction
	};

typedef btAlignedObjectArray<SContact> tSContactArray;
typedef btAlignedObjectArray<RContact> tRContactArray;
typedef btAlignedObjectArray<NetCollisionObject*> ObjArray;

public:
        Net_Wraper(net_t net, double P, double delta, double margin = 0.1);

        double getMargin() { return m_mrg; }
        void setMargin(double margin) { m_mrg = margin; }
        double computeFreeNexts();
        void addObjs2World(btCollisionWorld* wrld);
        void updateCollisionInfo();
        static void solveSContacts(Net_Wraper* body);
        static void solveRContacts(Net_Wraper* body);
        void setDelta(double delta) { m_delta = delta; }
        void setElasticModel(int type);
        void setP(double P) { m_P = P; }
private:
        void splitNet2ColObjs(const std::vector<int>& axisSeq, int depth);

private:
    ObjArray m_objs;
    tSContactArray m_scontacts;        // Soft contacts
    tRContactArray m_rcontacts;
	net_t m_net;
	double m_mrg = 0.1;
	double m_P;
	double m_delta;
	int m_eModel = 0;

	friend class World;
};

using NetObject = Net_Wraper::NetCollisionObject;

#endif // NET_WRAPER_H
