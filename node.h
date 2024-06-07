#ifndef node_h
#define node_h

/** \file node.h
\brief header to define struct Node
*/

//#include <memory>
#include <iostream>
#include <eigen3/Eigen/Dense>
#include "config.h"

/**
 \namespace Nodes
 to grab altogether some dedicated functions and data of the nodes
 */
namespace Nodes
    {
     /** space dimension */
    const int DIM = 3;
    
    /** size of the array of dataNode */
    const int NB_DATANODE = 2;

    /** convenient enum to avoid direct index manipulation in the array of struct dataNode */
    enum step
        {
        CURRENT = 0,
        NEXT = 1
        };

    /** convenient enum to specify what coordinate in calculations */
    enum index
        {
        IDX_UNDEF = -1,
        IDX_X = 0,
        IDX_Y = 1,
        IDX_Z = 2
        };

/** \return \f$ x^2 \f$ */
inline double sq(const double x) { return x * x; }

/** \struct dataNode
contains the vector fields u,v and the two associated scalar potentials
*/

struct dataNode
    {
    Eigen::Vector3d u;  /**< magnetization */
    Eigen::Vector3d v;  /**< magnetization speed */
    double phi;         /**< scalar potential */
    double phiv;        /**< scalar potential of velocity */
    };

/** \struct Node
Node is containing physical point of coordinates \f$ p = (x,y,z) \f$, magnetization value at \f$
m(p,t) \f$. Many other values for the computation of the scalar potential \f$ \phi \f$
*/
struct Node
    {
    Eigen::Vector3d p;  /**< Physical position p=(x,y,z)  of the node */
    Eigen::Vector3d ep; /**< local vector basis : \f$ e_p = \vec{rand} \times u0 \f$ , then normalized */
    Eigen::Vector3d eq; /**< local vector basis : \f$ e_q = u0 \times e_p \f$ , then normalized */
    
    /** datas associated to position p
    step CURRENT (0) : start of the time step
    step NEXT (1) : after the current time step 
    */ 
    dataNode d[NB_DATANODE];

    /** setter for the local basis vector */
    inline void setBasis(const double r)
        {
        // Choose for an initial ep the direction, among (X, Y, Z), which is further away from u0.
        Eigen::Index minIdx;
        d[CURRENT].u.cwiseAbs().minCoeff(&minIdx);
        #if EIGEN_VERSION_AT_LEAST(3,4,0)
            ep.setUnit(minIdx);
        #else
            if (minIdx == IDX_X) { ep = Eigen::Vector3d::UnitX(); }
            else if (minIdx == IDX_Y) { ep = Eigen::Vector3d::UnitY(); }
            else { ep = Eigen::Vector3d::UnitZ(); }
        #endif
        // Gram-Schmidt orthonormalization of (u0, ep).
        ep -= ep.dot(d[CURRENT].u) * d[CURRENT].u;
        ep.normalize();

        // Complete the basis with a vector product.
        eq = d[CURRENT].u.cross(ep);

        // Rotate (ep, eq) by the random angle.
        Eigen::Vector3d new_ep = cos(r) * ep - sin(r) * eq;
        eq = sin(r) * ep + cos(r) * eq;
        ep = new_ep;

        // The basis (u0, ep, eq) should already be orthonormal. An extra orthonormalization could
        // reduce the rounding errors.
        if (PARANOID_ORTHONORMALIZATION)
            {
            // Modified Gram-Schmidt orthonormalization.
            ep -= ep.dot(d[CURRENT].u) * d[CURRENT].u;
            ep.normalize();
            eq -= eq.dot(d[CURRENT].u) * d[CURRENT].u;
            eq -= eq.dot(ep) * ep;
            eq.normalize();
            }
        }

    /**
    preparation of the quantities u0,v0,phi0,phiv0 for incomming time-step
    */
    inline void evolution(void) { d[step::CURRENT] = d[step::NEXT]; }

    /**
    integration of the evolution of the magnetization for time step dt
    in a base composed of u0,ep,eq = u0*ep we have
    \f$ v = v_p e_p + v_q e_q \f$
    and new magnetization value is : \f$ u = u_0 + v dt \f$ after normalization
    */

    inline void make_evol(const double vp /**< [in] */, const double vq /**< [in] */,
                          const double dt /**< [in] */)
        {
        d[NEXT].v = vp * ep + vq * eq;
        d[NEXT].u = d[CURRENT].u + dt * d[NEXT].v;
        d[NEXT].u.normalize();
        }

    /** speed projection along ep */
    inline double proj_ep(void) const { return d[NEXT].v.dot(ep); }

    /** speed projection along eq */
    inline double proj_eq(void) const { return d[NEXT].v.dot(eq); }

    /** getter for u at step k */
    inline const Eigen::Vector3d get_u(step k /**< [in] */) const
    { return d[k].u; }

    /** getter for v at step k */
    inline const Eigen::Vector3d get_v(step k /**< [in] */) const 
    { return d[k].v; }

    };  // end struct node

/** getter for magnetizations */
template <step K>
Eigen::Vector3d get_u(Node const &n /**< [in] */) { return n.d[K].u; }

/** getter for speed magnetizations */
template <step K>
Eigen::Vector3d get_v(Node const &n /**< [in] */) { return n.d[K].v; }

/** getter for phis */
template <step K>
double get_phi(Node const &n /**< [in] */) { return n.d[K].phi; }

/** getter for phivs */
template <step K>
double get_phiv(Node const &n /**< [in] */) { return n.d[K].phiv; }

/** getter for u component */
inline double get_u_comp(Node const &n /**< [in] */, index idx /**< [in] */)
    { return n.d[NEXT].u(idx); }

/** getter for v component */
inline double get_v_comp(Node const &n /**< [in] */, index idx /**< [in] */)
    { return n.d[NEXT].v(idx); }

/** setter for phi */
inline void set_phi(Node &n, double val) { n.d[NEXT].phi = val; }

/** setter for phi_v */
inline void set_phiv(Node &n, double val) { n.d[NEXT].phiv = val; }

    }  // end namespace Nodes

#endif /* node_h */
