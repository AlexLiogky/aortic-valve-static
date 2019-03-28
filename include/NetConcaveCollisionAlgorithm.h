#ifndef NETCONCAVECOLLISIONALGORITHM_H
#define NETCONCAVECOLLISIONALGORITHM_H

#include "BulletCollision/BroadphaseCollision/btCollisionAlgorithm.h"
#include "BulletCollision/BroadphaseCollision/btDispatcher.h"
#include "BulletCollision/CollisionDispatch/btCollisionObjectWrapper.h"
#include "BulletCollision/CollisionShapes/btCollisionShape.h"
#include "BulletCollision/CollisionShapes/btConcaveShape.h"
#include "BulletCollision/CollisionDispatch/btCollisionObject.h"
#include "BulletCollision/BroadphaseCollision/btBroadphaseInterface.h"
#include "BulletCollision/CollisionShapes/btTriangleCallback.h"
#include "BulletCollision/NarrowPhaseCollision/btPersistentManifold.h"
class btDispatcher;
#include "BulletCollision/BroadphaseCollision/btBroadphaseProxy.h"
#include "BulletCollision/CollisionDispatch/btCollisionCreateFunc.h"
#include "Net_Wraper.h"
class btCollisionShape;

#include "LinearMath/btHashMap.h"

#include "BulletCollision/BroadphaseCollision/btQuantizedBvh.h"


struct btTriIndex
{
	int m_PartIdTriangleIndex;
	class btCollisionShape* m_childShape;

	btTriIndex(int partId, int triangleIndex, btCollisionShape* shape)
	{
		m_PartIdTriangleIndex = (partId << (31 - MAX_NUM_PARTS_IN_BITS)) | triangleIndex;
		m_childShape = shape;
	}

	int getTriangleIndex() const
	{
		// Get only the lower bits where the triangle index is stored
		unsigned int x = 0;
		unsigned int y = (~(x & 0)) << (31 - MAX_NUM_PARTS_IN_BITS);
		return (m_PartIdTriangleIndex & ~(y));
	}
	int getPartId() const
	{
		// Get only the highest bits where the part index is stored
		return (m_PartIdTriangleIndex >> (31 - MAX_NUM_PARTS_IN_BITS));
	}
	int getUid() const
	{
		return m_PartIdTriangleIndex;
	}
};

///For each triangle in the concave mesh that overlaps with the AABB of a soft body (m_softBody), processTriangle is called.
class NetTriangleCallback : public btTriangleCallback
{
	Net_Wraper* m_softBody;
	const btCollisionObject* m_triBody;

	btVector3 m_aabbMin;
	btVector3 m_aabbMax;

	btManifoldResult* m_resultOut;

	btDispatcher* m_dispatcher;
	const btDispatcherInfo* m_dispatchInfoPtr;
	btScalar m_collisionMarginTriangle;

	btHashMap<btHashKey<btTriIndex>, btTriIndex> m_shapeCache;

public:
	int m_triangleCount;

	NetTriangleCallback(btDispatcher* dispatcher, const btCollisionObjectWrapper* body0Wrap, const btCollisionObjectWrapper* body1Wrap, bool isSwapped);

	void setTimeStepAndCounters(btScalar collisionMarginTriangle, const btCollisionObjectWrapper* triObjWrap, const btDispatcherInfo& dispatchInfo, btManifoldResult* resultOut);

	virtual ~NetTriangleCallback() { clearCache(); }

	virtual void processTriangle(btVector3* triangle, int partId, int triangleIndex);

	void clearCache();

	SIMD_FORCE_INLINE const btVector3& getAabbMin() const
	{
		return m_aabbMin;
	}
	SIMD_FORCE_INLINE const btVector3& getAabbMax() const
	{
		return m_aabbMax;
	}
};



class NetConcaveCollisionAlgorithm : public btCollisionAlgorithm
{
	bool m_isSwapped;

	NetTriangleCallback m_NetTriangleCallback;

public:
	NetConcaveCollisionAlgorithm(const btCollisionAlgorithmConstructionInfo& ci, const btCollisionObjectWrapper* body0Wrap, const btCollisionObjectWrapper* body1Wrap, bool isSwapped)
	: btCollisionAlgorithm(ci),
	  m_isSwapped(isSwapped),
	  m_NetTriangleCallback(ci.m_dispatcher1, body0Wrap, body1Wrap, isSwapped)
    {}

	virtual ~NetConcaveCollisionAlgorithm() {};

	virtual void processCollision(const btCollisionObjectWrapper* body0Wrap, const btCollisionObjectWrapper* body1Wrap, const btDispatcherInfo& dispatchInfo, btManifoldResult* resultOut);

	btScalar calculateTimeOfImpact(btCollisionObject* body0, btCollisionObject* body1, const btDispatcherInfo& dispatchInfo, btManifoldResult* resultOut)
	{
        return 1.0f; //do Nothing
	}

	virtual void getAllContactManifolds(btManifoldArray& manifoldArray)
	{
		//we don't add any manifolds
	}

	void clearCache()
	{
        m_NetTriangleCallback.clearCache();
    }

	struct CreateFunc : public btCollisionAlgorithmCreateFunc
	{
		virtual btCollisionAlgorithm* CreateCollisionAlgorithm(btCollisionAlgorithmConstructionInfo& ci, const btCollisionObjectWrapper* body0Wrap, const btCollisionObjectWrapper* body1Wrap)
		{
			void* mem = btAlignedAlloc(static_cast<size_t>(sizeof(NetConcaveCollisionAlgorithm)), 16);//ci.m_dispatcher1->allocateCollisionAlgorithm(sizeof(NetConcaveCollisionAlgorithm));
			return new (mem) NetConcaveCollisionAlgorithm(ci, body0Wrap, body1Wrap, false);
		}
	};

	struct SwappedCreateFunc : public btCollisionAlgorithmCreateFunc
	{
		virtual btCollisionAlgorithm* CreateCollisionAlgorithm(btCollisionAlgorithmConstructionInfo& ci, const btCollisionObjectWrapper* body0Wrap, const btCollisionObjectWrapper* body1Wrap)
		{
			void* mem = btAlignedAlloc(static_cast<size_t>(sizeof(NetConcaveCollisionAlgorithm)), 16);
			//ci.m_dispatcher1->allocateCollisionAlgorithm(sizeof(NetConcaveCollisionAlgorithm));
			return new (mem) NetConcaveCollisionAlgorithm(ci, body0Wrap, body1Wrap, true);
		}
	};
};

#endif // NETCONCAVECOLLISIONALGORITHM_H
