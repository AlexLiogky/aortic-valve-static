#include "NetConcaveCollisionAlgorithm.h"
#include <iostream>

void NetConcaveCollisionAlgorithm::processCollision(const btCollisionObjectWrapper* body0Wrap, const btCollisionObjectWrapper* body1Wrap, const btDispatcherInfo& dispatchInfo, btManifoldResult* resultOut)
{
	//btCollisionObject* convexBody = m_isSwapped ? body1 : body0;
	const btCollisionObjectWrapper* triBody = m_isSwapped ? body0Wrap : body1Wrap;

	if (triBody->getCollisionShape()->isConcave())
	{
		const btCollisionObject* triOb = triBody->getCollisionObject();
		const btConcaveShape* concaveShape = static_cast<const btConcaveShape*>(triOb->getCollisionShape());

		//	if (convexBody->getCollisionShape()->isConvex())
		{
			btScalar collisionMarginTriangle = concaveShape->getMargin();

			//			resultOut->setPersistentManifold(m_btSoftBodyTriangleCallback.m_manifoldPtr);
			m_NetTriangleCallback.setTimeStepAndCounters(collisionMarginTriangle, triBody, dispatchInfo, resultOut);

			concaveShape->processAllTriangles(&m_NetTriangleCallback, m_NetTriangleCallback.getAabbMin(), m_NetTriangleCallback.getAabbMax());

			//	resultOut->refreshContactPoints();
		}
	}
}

NetTriangleCallback::NetTriangleCallback(btDispatcher* dispatcher, const btCollisionObjectWrapper* body0Wrap, const btCollisionObjectWrapper* body1Wrap, bool isSwapped) : m_dispatcher(dispatcher),
																																														 m_dispatchInfoPtr(0)
{
	m_softBody = (isSwapped ? (NetObject*)body1Wrap->getCollisionObject() : (NetObject*)body0Wrap->getCollisionObject());
	m_triBody = isSwapped ? body0Wrap->getCollisionObject() : body1Wrap->getCollisionObject();

	clearCache();
}

void NetTriangleCallback::clearCache()
{
	for (int i = 0; i < m_shapeCache.size(); i++)
	{
		btTriIndex* tmp = m_shapeCache.getAtIndex(i);
		btAssert(tmp);
		btAssert(tmp->m_childShape);
		//m_softBody->m_sparsesdf.RemoveReferences(tmp->m_childShape);  //necessary?
		delete tmp->m_childShape;
	}
	m_shapeCache.clear();
}

void NetTriangleCallback::processTriangle(btVector3* triangle, int partId, int triangleIndex)
{
	//std::cout << "I'm in processTriangle\n";
	/*btCollisionAlgorithmConstructionInfo ci;
	ci.m_dispatcher1 = m_dispatcher;

	btTriIndex triIndex(partId, triangleIndex, 0);
	btHashKey<btTriIndex> triKey(triIndex.getUid());

	btTriIndex* shapeIndex = m_shapeCache[triKey];*/
	m_softBody->CollisionHandler(m_triBody, triangle);
}

void NetTriangleCallback::setTimeStepAndCounters(btScalar collisionMarginTriangle, const btCollisionObjectWrapper* triBodyWrap, const btDispatcherInfo& dispatchInfo, btManifoldResult* resultOut)
{
	m_dispatchInfoPtr = &dispatchInfo;
	m_collisionMarginTriangle = collisionMarginTriangle + btScalar(0.06);
	m_resultOut = resultOut;

	btVector3 aabbWorldSpaceMin, aabbWorldSpaceMax;
	m_softBody->getAabb(aabbWorldSpaceMin, aabbWorldSpaceMax);
	btVector3 halfExtents = (aabbWorldSpaceMax - aabbWorldSpaceMin) * btScalar(0.5);
	btVector3 softBodyCenter = (aabbWorldSpaceMax + aabbWorldSpaceMin) * btScalar(0.5);

	btTransform softTransform;
	softTransform.setIdentity();
	softTransform.setOrigin(softBodyCenter);

	btTransform convexInTriangleSpace;
	convexInTriangleSpace = triBodyWrap->getWorldTransform().inverse() * softTransform;
	btTransformAabb(halfExtents, m_collisionMarginTriangle, convexInTriangleSpace, m_aabbMin, m_aabbMax);
}
