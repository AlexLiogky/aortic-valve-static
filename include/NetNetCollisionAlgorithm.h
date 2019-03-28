#ifndef NETNETCOLLISIONALGORITHM_H
#define NETNETCOLLISIONALGORITHM_H

#include "BulletCollision/BroadphaseCollision/btCollisionAlgorithm.h"
#include "BulletCollision/CollisionDispatch/btCollisionCreateFunc.h"
#include "BulletCollision/CollisionDispatch/btCollisionObjectWrapper.h"
#include "BulletCollision/BroadphaseCollision/btDispatcher.h"
#include "Net_Wraper.h"

class btPersistentManifold;

class NetNetCollisionAlgorithm : public btCollisionAlgorithm
{
    bool m_ownManifold;
	btPersistentManifold* m_manifoldPtr;

public:
	NetNetCollisionAlgorithm(const btCollisionAlgorithmConstructionInfo& ci)
		: btCollisionAlgorithm(ci) {}

	virtual void processCollision(const btCollisionObjectWrapper* body0Wrap, const btCollisionObjectWrapper* body1Wrap, const btDispatcherInfo& /*dispatchInfo*/, btManifoldResult* /*resultOut*/)
	{
        Net_Wraper* soft0 = (Net_Wraper*)body0Wrap->getCollisionObject();
        Net_Wraper* soft1 = (Net_Wraper*)body1Wrap->getCollisionObject();
        soft0->CollisionHandler(soft1);
    }

	virtual btScalar calculateTimeOfImpact(btCollisionObject* body0, btCollisionObject* body1, const btDispatcherInfo& dispatchInfo, btManifoldResult* resultOut)
	{
        return 1.f; //do Nothing
    }

	virtual void getAllContactManifolds(btManifoldArray& manifoldArray)
	{
		if (m_manifoldPtr && m_ownManifold)
			manifoldArray.push_back(m_manifoldPtr);
	}

	NetNetCollisionAlgorithm(btPersistentManifold* mf, const btCollisionAlgorithmConstructionInfo& ci, const btCollisionObjectWrapper* body0Wrap, const btCollisionObjectWrapper* body1Wrap)
        :NetNetCollisionAlgorithm(ci) {}

	virtual ~NetNetCollisionAlgorithm() {}

	struct CreateFunc : public btCollisionAlgorithmCreateFunc
	{
		virtual btCollisionAlgorithm* CreateCollisionAlgorithm(btCollisionAlgorithmConstructionInfo& ci, const btCollisionObjectWrapper* body0Wrap, const btCollisionObjectWrapper* body1Wrap)
		{
			int bbsize = sizeof(NetNetCollisionAlgorithm);
			void* ptr = ci.m_dispatcher1->allocateCollisionAlgorithm(bbsize);
			return new (ptr) NetNetCollisionAlgorithm(0, ci, body0Wrap, body1Wrap);
		}
	};
};

#endif // NETNETCOLLISIONALGORITHM_H
