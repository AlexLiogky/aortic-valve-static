#ifndef STATICCOLLISIONALGORITHM_H
#define STATICCOLLISIONALGORITHM_H

#include "BulletCollision/CollisionShapes/btTriangleCallback.h"
#include "BulletCollision/CollisionDispatch/btCollisionObject.h"
#include "Net_Wraper.h"

///For each triangle in the concave mesh that overlaps with the AABB of a soft body (m_softBody), processTriangle is called.
class Net_WraperTriangleCallback : public btTriangleCallback
{
	NetObject* m_softBody;
	const btCollisionObject* m_triBody;

	btVector3 m_aabbMin;
	btVector3 m_aabbMax;

	//btManifoldResult* m_resultOut;

	double m_collisionMarginTriangle;

	//btHashMap<btHashKey<btTriIndex>, btTriIndex> m_shapeCache;

public:
	//int m_triangleCount;

	//	btPersistentManifold*	m_manifoldPtr;

	Net_WraperTriangleCallback(NetObject* body0, const btCollisionObject* body1):
        m_softBody{body0},
        m_triBody{body1}
    {

    }

	//virtual ~btSoftBodyTriangleCallback();

	virtual void processTriangle(btVector3* triangle, int partId, int triangleIndex);

	//void clearCache();

	SIMD_FORCE_INLINE const btVector3& getAabbMin() const
	{
		return m_aabbMin;
	}
	SIMD_FORCE_INLINE const btVector3& getAabbMax() const
	{
		return m_aabbMax;
	}

	void setSizes(double collisionMarginTriangle);
};


#endif // STATICCOLLISIONALGORITHM_H
