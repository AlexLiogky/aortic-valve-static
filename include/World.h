#ifndef WORLD_H
#define WORLD_H

#include "computation.h"
#include "nets.h"
#include <vector>
#include "Net_Wraper.h"
#include "BulletCollision/CollisionShapes/btCollisionShape.h"
#include "BulletCollision/CollisionShapes/btTriangleIndexVertexArray.h"
#include "BulletCollision/CollisionShapes/btBvhTriangleMeshShape.h"
#include "BulletCollision/CollisionDispatch/btCollisionObject.h"
#include "StaticCollisionAlgorithm.h"
#include "BulletCollision/CollisionDispatch/btDefaultCollisionConfiguration.h"
#include "BulletCollision/CollisionDispatch/btCollisionWorld.h"
#include "BulletCollision/BroadphaseCollision/btAxisSweep3.h"
#include "NetNetCollisionAlgorithm.h"
#include "NetConcaveCollisionAlgorithm.h"

using namespace std;

class World
{
public:
        nets_t& getDynamicNets() { return m_dynamic_nets; }
        nets_t& getStaticNets() { return m_static_nets; }
        World(nets_t& dynamic_nets, nets_t& static_nets, wrld_cnd_t& cond, solver_t& solver_data, double drop_thr, double max_shft);
        void compute_nets_time(long double compute_time = 1000.0/60, int max_its = 10000);

private:
    void set_initial_solving_params();
    void convert_net_to_btTriangleMesh(net_t& st);
    //void findSoftCollisions();
    //void findRigidCollisions();
    //void solveStaticCollision(btCollisionObject* b0, Net_Wraper* b1);
    void addTestBvh();
    void registerCollisionWorld();
protected:
    double predictMotion();
    void findCollisions();
    void solveCollisions();

private:
    nets_t m_dynamic_nets;
    nets_t m_static_nets;
    nets_t m_union_nets;
    std::vector<Net_Wraper*> m_collision;
    std::vector<btCollisionObject*> m_static_col;
    btCollisionWorld* m_collisionWorld;
    wrld_cnd_t m_conditions;
    solver_t m_solver_data;
    statistic_t m_statistic;
	double m_allow_shift = 1;
	double m_max_shift;

};

class MyCollisionConfiguration : public btDefaultCollisionConfiguration
{
    btCollisionAlgorithmCreateFunc* m_softSoftCreateFunc;
    btCollisionAlgorithmCreateFunc* m_softBodyConcaveCreateFunc;
	btCollisionAlgorithmCreateFunc* m_swappedSoftBodyConcaveCreateFunc;
public:
    MyCollisionConfiguration(const btDefaultCollisionConstructionInfo& constructionInfo = btDefaultCollisionConstructionInfo());
    virtual ~MyCollisionConfiguration();
    ///creation of soft-soft and soft-rigid, and otherwise fallback to base class implementation
	virtual btCollisionAlgorithmCreateFunc* getCollisionAlgorithmCreateFunc(int proxyType0, int proxyType1) override;
};


#endif // WORLD_H
