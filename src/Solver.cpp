#include "Solver.h"

Solver::Solver(nets_t& dynamic_nets, nets_t& static_nets, wrld_cnd_t& cond, solver_t& solver_data, double drop_thr, double max_shft):
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

    set_initial_solving_params();
    for (unsigned int i = 0; i < m_static_nets.count; ++i)
        convert_net_to_btTriangleMesh(m_static_nets.nets[i]);
}


void Solver::convert_net_to_btTriangleMesh(const net_t& st)
{
    const int totalTriangles = st.elems.count;
    const int totalVerts = st.vrtx.count;
    btVector3* vertices = new btVector3[totalVerts];
    for (int i = 0; i < totalVerts; ++i)
        vertices[i] = p2Bt(st.vrtx.nodes[i]->coord);
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
                                                                                    vertStride);
    bool useQuantizedAabbCompression = true;
    btCollisionShape* meshShape = new btBvhTriangleMeshShape(indexVertexArrays, useQuantizedAabbCompression);
    meshShape->setMargin(0.2);
    btCollisionObject* newOb = new btCollisionObject();
    btTransform tr;
	tr.setIdentity();
    newOb->setWorldTransform(tr);
	newOb->setInterpolationWorldTransform(tr);
    newOb->setCollisionShape(meshShape);
    m_static_col.push_back(newOb);
}

void Solver::set_initial_solving_params()
{
    predictMotion();
    findSoftCollisions();
    solveCollisions();
    update_statistic_t(m_dynamic_nets, &m_statistic);

    m_statistic.init_div = statistic_t_get_full_mid_diviation(m_statistic);
    statistic_t_reset(&m_statistic);
}

double Solver::predictMotion(){
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

void Solver::findSoftCollisions(){
    for (unsigned int i = 1; i < m_collision.size(); ++i)
        for (unsigned int j = 0; j < i; ++j)
        {
            m_collision[i]->CollisionHandler(m_collision[j]);
        }
}

void Solver::findRigidCollisions(){
    for (unsigned int i = 1; i < m_collision.size(); ++i)
        for (unsigned int j = 0; j < m_static_col.size(); ++j)
        {
            solveStaticCollision(m_static_col[j], m_collision[i]);
        }
}

void Solver::findCollisions()
{
    findSoftCollisions();
    findRigidCollisions();
}

void Solver::solveCollisions(){
    for (auto& body: m_collision)
    {
        Net_Wraper::solveSContacts(body);
        Net_Wraper::solveRContacts(body);
    }
}

void Solver::compute_nets_time(long double compute_time, int max_its)
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

void Solver::solveStaticCollision(btCollisionObject* b0, Net_Wraper* b1)
{
    Net_WraperTriangleCallback callback(b1, b0);
    const btConcaveShape* concaveShape = static_cast<const btConcaveShape*>(b0->getCollisionShape());
    btScalar collisionMarginTriangle = concaveShape->getMargin();
    callback.setSizes(collisionMarginTriangle);
    concaveShape->processAllTriangles(&callback, callback.getAabbMin(), callback.getAabbMax());
}


