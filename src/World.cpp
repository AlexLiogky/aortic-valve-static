#include "World.h"
#include <iostream>
#include <algorithm>
#include "bound-box.h"


World::World(nets_t& dynamic_nets, nets_t& static_nets, wrld_cnd_t& cond, solver_t& solver_data, double drop_thr, double max_shft):
    m_dynamic_nets{dynamic_nets},
    m_static_nets{static_nets},
    m_union_nets{create_union_net(dynamic_nets, static_nets)},
    m_conditions{cond},
    m_solver_data{solver_data},
    m_statistic{statistical_data_construct(dynamic_nets)},
    m_allow_shift{drop_thr},
    m_max_shift{max_shft}
{
    m_collision.reserve(m_union_nets.count);
    for (unsigned int i = 0; i < m_dynamic_nets.count; ++i)
        m_collision.push_back(new Net_Wraper(m_dynamic_nets.nets[i], m_conditions.P, m_solver_data.delta));

    point_t_dump(Bt2p(m_collision[0]->m_bounds[0]));
    point_t_dump(Bt2p(m_collision[0]->m_bounds[1]));


    for (unsigned int i = 0; i < m_static_nets.count; ++i)
        convert_net_to_btTriangleMesh(m_static_nets.nets[i]);

    //addTestBvh();

    registerCollisionWorld();
    set_initial_solving_params();
}

void World::registerCollisionWorld()
{
    btVector3 worldAabbMin(-100, -100, -100);
	btVector3 worldAabbMax(100, 100, 100);
    const int maxProxies = 1000;

	btBroadphaseInterface* broadphase = new btAxisSweep3(worldAabbMin, worldAabbMax, maxProxies);
    btCollisionConfiguration* collisionConfig = new MyCollisionConfiguration(btDefaultCollisionConstructionInfo());
    btDispatcher* dispatcher = new btCollisionDispatcher(collisionConfig);
    m_collisionWorld = new btCollisionWorld(dispatcher, broadphase, collisionConfig);

    for (auto& i: m_static_col)
        m_collisionWorld->addCollisionObject(i, btBroadphaseProxy::StaticFilter, btBroadphaseProxy::AllFilter ^ btBroadphaseProxy::StaticFilter);
    for (auto& i: m_collision)
        m_collisionWorld->addCollisionObject(i, btBroadphaseProxy::DefaultFilter, btBroadphaseProxy::AllFilter);
}

bool inline point_t_in_Aabb(const point_t& p, const double brds[3][2])
{
    bool res = true;
    for (int i = 0; i < 3; ++i)
        res &= (p.coord[i] >= brds[i][0]) && (p.coord[i] <= brds[i][1]);
    return res;
}

btTriangleIndexVertexArray* split_net_to_Aabb(const net_t& net, const double brds[3][2])
{
    int* vrt_in_Aabb = new int[net.vrtx.count];
    std::fill<int*, int>(vrt_in_Aabb, vrt_in_Aabb + net.vrtx.count, -1);
    int totalTriangles = 0;
    for (int i = 0, nt = net.elems.count; i < nt; ++i)
    {
        node_t** e = net.elems.elems[i]->vrts;
        int check = 0;
        for (int j = 0; j < 3; ++j)
        {
            if (point_t_in_Aabb(e[j]->coord, brds))
            {
                vrt_in_Aabb[e[j]->id] = 1;
                ++check;
            }
        }
        totalTriangles += (check == 3);
    }
    int totalVerts = 0;
    for (int i = 0, nv = net.vrtx.count; i < nv; ++i)
        totalVerts += (vrt_in_Aabb[i] == 1);
    btVector3* vertices = new btVector3[totalVerts];
    for (int i = 0, j = 0, nv = net.vrtx.count; i < nv; ++i)
        if (vrt_in_Aabb[i] == 1)
        {
            vrt_in_Aabb[i] = j;
            vertices[j++] = p2Bt(net.vrtx.nodes[i]->coord);
        }

    int* triag_indeces = new int[totalTriangles * 3];
    for (int i = 0, k = 0, nt = net.elems.count; i < nt; ++i)
    {
        elem_t& e = *net.elems.elems[i];
        bool in = true;
        for (int j = 0; j < 3; ++j)
             in &= (vrt_in_Aabb[e.vrts[j]->id] != -1);
        if (!in) continue;
        if (e.coef >= 0)
            for (int j = 0; j < 3; ++j)
               triag_indeces[3 * k + j] = vrt_in_Aabb[e.vrts[j]->id];
        else
            for (int j = 0; j < 3; ++j)
               triag_indeces[3 * k + 2 - j] = vrt_in_Aabb[e.vrts[j]->id];
        ++k;
    }

    delete[] vrt_in_Aabb;

    const int vertStride = sizeof(btVector3);
    const int indexStride = 3 * sizeof(int);
    return new btTriangleIndexVertexArray(  totalTriangles, triag_indeces, indexStride,
                                            totalVerts, (btScalar*)&vertices[0].x(),
                                            vertStride                                  );
}

void World::convert_net_to_btTriangleMesh(net_t& st)
{
   /* const int totalTriangles = st.elems.count;
    const int totalVerts = st.vrtx.count;
    btVector3* vertices = new btVector3[totalVerts];
    for (int i = 0; i < totalVerts; ++i){
        vertices[i] = p2Bt(st.vrtx.nodes[i]->coord);
        st.vrtx.nodes[i]->id = i;
    }
    int* triag_indeces = new int[totalTriangles * 3];
    for (int i = 0; i < totalTriangles; ++i)
    {
        elem_t& e = *st.elems.elems[i];
        if (e.coef >= 0)
            for (int j = 0; j < 3; ++j)
               triag_indeces[3 * i + j] = e.vrts[j]->id;
        else
            for (int j = 0; j < 3; ++j)
               triag_indeces[3 * i + 2 - j] = e.vrts[j]->id;
    }
    const int vertStride = sizeof(btVector3);
    const int indexStride = 3 * sizeof(int);
    btTriangleIndexVertexArray* indexVertexArrays = new btTriangleIndexVertexArray( totalTriangles,
                                                                                    triag_indeces,
                                                                                    indexStride,
                                                                                    totalVerts,
                                                                                    (btScalar*)&vertices[0].x(),
                                                                                    vertStride);*/
    double brds[3][2];
    nets_t_fill_borders(m_dynamic_nets, brds);
	double ext_coef = 1.1;
	extend_borders(ext_coef, brds);
	btTriangleIndexVertexArray* indexVertexArrays = split_net_to_Aabb(st, brds);
    bool useQuantizedAabbCompression = true;
    btBvhTriangleMeshShape* meshShape = new btBvhTriangleMeshShape(indexVertexArrays, useQuantizedAabbCompression);
    meshShape->setMargin(0.5);//16);
    meshShape->buildOptimizedBvh();
    btCollisionObject* newOb = new btCollisionObject();
    btTransform tr;
	tr.setIdentity();
    newOb->setWorldTransform(tr);
	newOb->setInterpolationWorldTransform(tr);
    newOb->setCollisionShape(meshShape);
    m_static_col.push_back(newOb);
}

void World::addTestBvh()
{
    const float TRIANGLE_SIZE = 1.0f;
    static float waveheight = 5.f;
    btCollisionShape* groundShape = 0;
	{
		int i;
		int j;

		const int NUM_VERTS_X = 300;
		const int NUM_VERTS_Y = 300;
		const int totalVerts = NUM_VERTS_X * NUM_VERTS_Y;
		const int totalTriangles = 2 * (NUM_VERTS_X - 1) * (NUM_VERTS_Y - 1);

		btVector3* gGroundVertices = new btVector3[totalVerts];
		int* gGroundIndices = new int[totalTriangles * 3];

		btScalar offset(0);

		for (i = 0; i < NUM_VERTS_X; i++)
		{
			for (j = 0; j < NUM_VERTS_Y; j++)
			{
				gGroundVertices[i + j * NUM_VERTS_X].setValue((i - NUM_VERTS_X * 0.5f) * TRIANGLE_SIZE,
															  0.f,
															  //waveheight * sinf((float)i/100) * cosf((float)j/100 + offset),
															  (j - NUM_VERTS_Y * 0.5f) * TRIANGLE_SIZE);
			}
		}

		int vertStride = sizeof(btVector3);
		int indexStride = 3 * sizeof(int);

		int index = 0;
		for (i = 0; i < NUM_VERTS_X - 1; i++)
		{
			for (int j = 0; j < NUM_VERTS_Y - 1; j++)
			{
				gGroundIndices[index++] = j * NUM_VERTS_X + i;
				gGroundIndices[index++] = (j + 1) * NUM_VERTS_X + i + 1;
				gGroundIndices[index++] = j * NUM_VERTS_X + i + 1;
				;

				gGroundIndices[index++] = j * NUM_VERTS_X + i;
				gGroundIndices[index++] = (j + 1) * NUM_VERTS_X + i;
				gGroundIndices[index++] = (j + 1) * NUM_VERTS_X + i + 1;
			}
		}

		btTriangleIndexVertexArray* indexVertexArrays = new btTriangleIndexVertexArray(totalTriangles,
																					   gGroundIndices,
																					   indexStride,
																					   totalVerts, (btScalar*)&gGroundVertices[0].x(), vertStride);

		bool useQuantizedAabbCompression = false;

		groundShape = new btBvhTriangleMeshShape(indexVertexArrays, useQuantizedAabbCompression);
		groundShape->setMargin(10);
	}
	btTransform tr;
	tr.setIdentity();
	tr.setOrigin(btVector3(-5, -20, 10));

	btCollisionObject* newOb = new btCollisionObject();
	newOb->setWorldTransform(tr);
	newOb->setInterpolationWorldTransform(tr);
	newOb->setCollisionShape(groundShape);
	m_static_col.push_back(newOb);
}

void World::set_initial_solving_params()
{
    predictMotion();
    findCollisions();
    solveCollisions();
    update_statistic_t(m_dynamic_nets, &m_statistic);

    m_statistic.init_div = statistic_t_get_full_mid_diviation(m_statistic);
    statistic_t_reset(&m_statistic);
}

double World::predictMotion(){
    double shift = 0;
    for (auto& net: m_collision){
        double max_net_shift = net->computeFreeNexts();
        if (shift < max_net_shift) shift = max_net_shift;
    }
    char ret = 0;
    if (shift > m_max_shift) {
		relaxation(m_union_nets, sqrt(m_allow_shift / shift));
		ret++;
	}
    for (auto& net: m_collision){
        net->updateCollisionInfo();
    }
    return ret;
}

void World::solveCollisions()
{
    for (auto& body: m_collision)
    {
        Net_Wraper::solveSContacts(body);
        Net_Wraper::solveRContacts(body);
    }
}

void World::findCollisions()
{
    m_collisionWorld->performDiscreteCollisionDetection();
}

void World::compute_nets_time(long double compute_time, int max_its)
{
    struct timeval start, end;
	gettimeofday(&start, NULL);
	int i = 0, crush = 0;
	for (i = 0; i < max_its; ++i){
        //collision_data_t_update(m_union_nets, &world->collision);
        crush += predictMotion();
        findCollisions();
        solveCollisions();

        update_statistic_t(m_dynamic_nets, &m_statistic);
        update_nets(m_dynamic_nets);

        gettimeofday(&end, NULL);
		if (compute_dif_ms(end, start) >= compute_time - 0.8) break;
	}
	double res = statistic_t_get_full_mid_diviation(m_statistic);
	double rms = res / m_statistic.cnt_nodes;
	//printf("RMS shift per iter of node = %e mm, relation = %lg, cr = %d / %d\n", rms , res / world->statistic.init_div, crush, i+1);
	statistic_t_reset(&m_statistic);
}

MyCollisionConfiguration::MyCollisionConfiguration(const btDefaultCollisionConstructionInfo& constructionInfo)
    : btDefaultCollisionConfiguration(constructionInfo)
{
    void* mem;
    mem = btAlignedAlloc(sizeof(NetNetCollisionAlgorithm::CreateFunc), 16);
	m_softSoftCreateFunc = new (mem) NetNetCollisionAlgorithm::CreateFunc;

	mem = btAlignedAlloc(sizeof(NetConcaveCollisionAlgorithm::CreateFunc), 16);
	m_softBodyConcaveCreateFunc = new (mem) NetConcaveCollisionAlgorithm::CreateFunc;

	mem = btAlignedAlloc(sizeof(NetConcaveCollisionAlgorithm::CreateFunc), 16);
	m_swappedSoftBodyConcaveCreateFunc = new (mem) NetConcaveCollisionAlgorithm::SwappedCreateFunc;
	m_swappedSoftBodyConcaveCreateFunc->m_swapped = true;
}

MyCollisionConfiguration::~MyCollisionConfiguration()
{
	m_softSoftCreateFunc->~btCollisionAlgorithmCreateFunc();
	btAlignedFree(m_softSoftCreateFunc);
    m_softBodyConcaveCreateFunc->~btCollisionAlgorithmCreateFunc();
	btAlignedFree(m_softBodyConcaveCreateFunc);
	m_swappedSoftBodyConcaveCreateFunc->~btCollisionAlgorithmCreateFunc();
	btAlignedFree(m_swappedSoftBodyConcaveCreateFunc);
}

btCollisionAlgorithmCreateFunc* MyCollisionConfiguration::getCollisionAlgorithmCreateFunc(int proxyType0, int proxyType1)
{
    ///softbody versus softbody
    if ((proxyType0 == SOFTBODY_SHAPE_PROXYTYPE) && (proxyType1 == SOFTBODY_SHAPE_PROXYTYPE))
	{
		return m_softSoftCreateFunc;
	}

	///softbody versus convex
	if (proxyType0 == SOFTBODY_SHAPE_PROXYTYPE && btBroadphaseProxy::isConcave(proxyType1))
	{
		return m_softBodyConcaveCreateFunc;
	}

	///convex versus soft body
	if (btBroadphaseProxy::isConcave(proxyType0) && proxyType1 == SOFTBODY_SHAPE_PROXYTYPE)
	{
		return m_swappedSoftBodyConcaveCreateFunc;
	}
    return btDefaultCollisionConfiguration::getCollisionAlgorithmCreateFunc(proxyType0, proxyType1);;
}

