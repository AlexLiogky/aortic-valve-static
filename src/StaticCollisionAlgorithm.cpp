#include "StaticCollisionAlgorithm.h"
#include "BulletCollision/CollisionShapes/btCollisionShape.h"
#include "BulletCollision/CollisionDispatch/btCollisionObject.h"


void Net_WraperTriangleCallback::setSizes(double collisionMarginTriangle)
{
	m_collisionMarginTriangle = collisionMarginTriangle + 0.06;

	btVector3 aabbWorldSpaceMin, aabbWorldSpaceMax;
	m_softBody->getAabb(aabbWorldSpaceMin, aabbWorldSpaceMax);
	btVector3 halfExtents = (aabbWorldSpaceMax - aabbWorldSpaceMin) * btScalar(0.5);
	btVector3 softBodyCenter = (aabbWorldSpaceMax + aabbWorldSpaceMin) * btScalar(0.5);
    btTransform softTransform;
	softTransform.setIdentity();
	softTransform.setOrigin(softBodyCenter);

	btTransformAabb(halfExtents, m_collisionMarginTriangle, softTransform, m_aabbMin, m_aabbMax);
}


void Net_WraperTriangleCallback::processTriangle(btVector3* triangle, int partId, int triangleIndex)
{
		m_softBody->CollisionHandler(m_triBody, triangle);
}
