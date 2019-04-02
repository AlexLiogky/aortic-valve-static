#ifndef NETCOLLISIONSHAPE_H
#define NETCOLLISIONSHAPE_H

#include "BulletCollision/CollisionShapes/btCollisionShape.h"
#include "BulletCollision/CollisionShapes/btConcaveShape.h"
#include "LinearMath/btQuickprof.h"
#include "LinearMath/btPolarDecomposition.h"
#include "BulletCollision/BroadphaseCollision/btBroadphaseInterface.h"
#include "BulletCollision/CollisionDispatch/btCollisionDispatcher.h"
#include "BulletCollision/CollisionShapes/btConvexInternalShape.h"
#include "BulletCollision/NarrowPhaseCollision/btGjkEpa2.h"

#include "Net_Wraper.h"

class NetCollisionShape : public btConcaveShape
{
    public:
    NetObject* m_body;

    NetCollisionShape(NetObject* backptr)
    {
		m_shapeType = SOFTBODY_SHAPE_PROXYTYPE;
		m_body = backptr;
	}
    virtual ~NetCollisionShape() {}

    void processAllTriangles(btTriangleCallback* /*callback*/, const btVector3& /*aabbMin*/, const btVector3& /*aabbMax*/) const
	{
		//not yet
		assert(0 && "processAllTriangles");
	}
	virtual void setLocalScaling(const btVector3& /*scaling*/)
	{
		///na
	}
	virtual const btVector3& getLocalScaling() const
	{
		static const btVector3 dummy(1, 1, 1);
		return dummy;
	}
	virtual void calculateLocalInertia(btScalar /*mass*/, btVector3& /*inertia*/) const
	{
		///not yet
		assert(0 && "calculateLocalInertia");
	}
	virtual const char* getName() const
	{
		return "NetBody";
	}

    virtual void getAabb(const btTransform& t, btVector3& aabbMin, btVector3& aabbMax) const
    {
        //m_body->getAabb(aabbMin, aabbMax);
        const btVector3 mins= m_body->m_bounds[0];
        const btVector3 maxs= m_body->m_bounds[1];
        const btVector3 crns[]={t*btVector3(mins.x(),mins.y(),mins.z()),
                         t*btVector3(maxs.x(),mins.y(),mins.z()),
                         t*btVector3(maxs.x(),maxs.y(),mins.z()),
                         t*btVector3(mins.x(),maxs.y(),mins.z()),
                         t*btVector3(mins.x(),mins.y(),maxs.z()),
                         t*btVector3(maxs.x(),mins.y(),maxs.z()),
                         t*btVector3(maxs.x(),maxs.y(),maxs.z()),
                         t*btVector3(mins.x(),maxs.y(),maxs.z())};
         aabbMin=aabbMax=crns[0];
         for(int i=1;i<8;++i)
         {
                 aabbMin.setMin(crns[i]);
                 aabbMax.setMax(crns[i]);
         }
    }


    protected:

    private:
};

#endif // NETCOLLISIONSHAPE_H
