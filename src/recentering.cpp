#include <functional>

#include "fem.h"

using namespace Nodes;

void Fem::direction(enum index idx_dir)
    {
    const int NPS = 1;

    ANNidxArray nnIdx = new ANNidx[NPS];
    if (!nnIdx) SYSTEM_ERROR;
    ANNdistArray dists = new ANNdist[NPS];
    if (!dists) SYSTEM_ERROR;
    ANNpoint queryPt = annAllocPt(DIM);

    Eigen::Vector3d qPt = msh.c - 0.5 * msh.l.cwiseProduct(Eigen::Vector3d::Unit(idx_dir));

    /* left */
    queryPt[0] = qPt.x();
    queryPt[1] = qPt.y();
    queryPt[2] = qPt.z();

    kdtree->annkSearch(queryPt, NPS, nnIdx, dists, 0.);
    int ns = nnIdx[0];
    double u2L = msh.getNode_u(ns)(idx_dir);

    /* right */
    qPt = msh.c + 0.5 * msh.l.cwiseProduct(Eigen::Vector3d::Unit(idx_dir));
    queryPt[0] = qPt.x();
    queryPt[1] = qPt.y();
    queryPt[2] = qPt.z();

    kdtree->annkSearch(queryPt, NPS, nnIdx, dists, 0.);
    ns = nnIdx[0];
    double u2R = msh.getNode_u(ns)(idx_dir);

    if (u2L * u2R > 0.)
        DW_dir = 0.0;
    else
        DW_dir = (u2L > 0 ? 1. : -1.); /* direction of the DW propagation +Oz or -Oz */
    delete[] dists;
    delete[] nnIdx;
    annDeallocPt(queryPt);
    }

bool Fem::recenter(double thres, char recentering_direction)
    {
    enum index idx_dir;

    switch (recentering_direction)
        {
        case 'X': idx_dir = IDX_X; break;
        case 'Y': idx_dir = IDX_Y; break;
        case 'Z': idx_dir = IDX_Z; break;
        default: idx_dir = IDX_Z; break;
        }

    thres = abs(thres);
    double m_i = msh.avg(Nodes::get_u_comp, idx_dir);
    if (fabs(m_i) < thres) return false;

    const int NPS = 1;
    int ns;

    ANNidxArray nnIdx = new ANNidx[NPS];
    if (!nnIdx) SYSTEM_ERROR;
    ANNdistArray dists = new ANNdist[NPS];
    if (!dists) SYSTEM_ERROR;
    ANNpoint queryPt = annAllocPt(DIM);

    Eigen::Vector3d p_dir = Eigen::Vector3d::Unit(idx_dir);  // unit vector
    Eigen::Vector3d qPt = msh.c - 0.5 * p_dir.cwiseProduct(msh.l);

    /* left frontier */
    queryPt[0] = qPt.x();
    queryPt[1] = qPt.y();
    queryPt[2] = qPt.z();

    kdtree->annkSearch(queryPt, NPS, nnIdx, dists, 0.);
    ns = nnIdx[0];
    double u2L = msh.getNode_u(ns)(idx_dir);

    // right frontier
    qPt = msh.c + 0.5 * p_dir.cwiseProduct(msh.l);
    queryPt[0] = qPt.x();
    queryPt[1] = qPt.y();
    queryPt[2] = qPt.z();

    kdtree->annkSearch(queryPt, NPS, nnIdx, dists, 0.);
    ns = nnIdx[0];
    double u2R = msh.getNode_u(ns)(idx_dir);

    if (u2L * u2R > 0)
        {
        std::cout << "Error No Domain Wall" << std::endl;
        SYSTEM_ERROR;
        }

    assert(u2L * u2R < 0);
    double D_i = m_i * 0.5 * msh.l(idx_dir) * u2L;  // shift with correct sign

    for (int i = 0; i < msh.getNbNodes(); i++)
        {
        Eigen::Vector3d p = msh.getNode_p(i) + D_i * p_dir;

        if (p(idx_dir) - msh.c(idx_dir) > +0.5 * msh.l(idx_dir))
            {
            p = msh.c + 0.5 * msh.l(idx_dir) * p_dir;
            }
        if (p(idx_dir) - msh.c(idx_dir) < -0.5 * msh.l(idx_dir))
            {
            p = msh.c - 0.5 * msh.l(idx_dir) * p_dir;
            }

        queryPt[0] = p.x();
        queryPt[1] = p.y();
        queryPt[2] = p.z();
        kdtree->annkSearch(queryPt, NPS, nnIdx, dists, 0.);

        int ns = nnIdx[0];
        msh.set_node_u0(i, msh.getNode_u(ns));
        msh.set_node_zero_v(i);
        }

    annDeallocPt(queryPt);
    delete[] nnIdx;
    delete[] dists;

    return true;
    }
