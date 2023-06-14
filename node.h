#ifndef node_h
#define node_h

/** \file node.h
\brief header to define struct Node
*/

#include <memory>

#include "config.h"
#include "pt3D.h"

/**
 \namespace Nodes
 to grab altogether some dedicated functions and data of the nodes
 */
namespace Nodes
    {
/**
 *convenient enum mainly to access some data values in calculations, used by some templates
 */
enum index
    {
    IDX_p = 0,
    IDX_u0 = 1,
    IDX_v0 = 2,
    IDX_u = 3,
    IDX_v = 4,
    IDX_phi0 = 5,
    IDX_phi = 6,
    IDX_phiv0 = 7,
    IDX_phiv = 8
    };

/** \struct Node
Node is containing physical point of coordinates \f$ p = (x,y,z) \f$, magnetization value at \f$
m(p,t) \f$. Many other values for the computation of the scalar potential \f$ \phi \f$
*/
struct Node
    {
    Pt::pt3D p;  /**< Physical position p=(x,y,z)  of the node */
    Pt::pt3D u0; /**< magnetization at the start of the current time step */
    Pt::pt3D v0; /**< magnetization speed at the start of the current time step */
    Pt::pt3D u;  /**< magnetization after the current time step */
    Pt::pt3D v;  /**< magnetization speed after the current time step */

    Pt::pt3D ep; /**< local vector basis : \f$ e_p = \vec{rand} \times u0 \f$ , then normalized */
    Pt::pt3D eq; /**< local vector basis : \f$ e_q = u0 \times e_p \f$ , then normalized */

    double phi0;  /**< scalar potential at the start of the current time step */
    double phi;   /**< scalar potential after the current time step */
    double phiv0; /**< scalar potential of velocity at the start of the current time step */
    double phiv;  /**< scalar potential of velocity after the current time step */

    /** setter for the local basis vector */
    inline void setBasis(const double r)
        {
        // Choose for an initial ep the direction, among (X, Y, Z), which is further away from u0.
        double abs_x = fabs(u0.x()), abs_y = fabs(u0.y()), abs_z = fabs(u0.z());
        if (abs_x < abs_y)
            {
            ep = abs_x < abs_z ? Pt::pt3D(Pt::IDX_X) : Pt::pt3D(Pt::IDX_Z);
            }
        else
            {
            ep = abs_y < abs_z ? Pt::pt3D(Pt::IDX_Y) : Pt::pt3D(Pt::IDX_Z);
            }

        // Gram-Schmidt orthonormalization of (u0, ep).
        ep -= pScal(ep, u0) * u0;
        ep.normalize();

        // Complete the basis with a vector product.
        eq = u0 * ep;

        // Rotate (ep, eq) by the random angle.
        Pt::pt3D new_ep = cos(r) * ep - sin(r) * eq;
        eq = sin(r) * ep + cos(r) * eq;
        ep = new_ep;

        // The basis (u0, ep, eq) should already be orthonormal. An extra orthonormalization could
        // reduce the rounding errors.
        if (PARANOID_ORTHONORMALIZATION)
            {
            // Modified Gram-Schmidt orthonormalization.
            ep -= pScal(ep, u0) * u0;
            ep.normalize();
            eq -= pScal(eq, u0) * u0;
            eq -= pScal(eq, ep) * ep;
            eq.normalize();
            }
        }

    /**
    preparation of the quantities u0,v0,phi0,phiv0 for incomming time-step
    */
    inline void evolution(void)
        {
        u0 = u;
        v0 = v;
        phi0 = phi;
        phiv0 = phiv;
        }

    /**
    integration of the evolution of the magnetization for time step dt
    in a base composed of u0,ep,eq = u0*ep we have
    \f$ v = v_p e_p + v_q e_q \f$
    and new magnetization value is : \f$ u = u_0 + v dt \f$ after normalization
    */

    inline void make_evol(const double vp /**< [in] */, const double vq /**< [in] */,
                          const double dt /**< [in] */)
        {
        v = vp * ep + vq * eq;
        u = u0 + dt * v;
        u.normalize();
        }

    };  // end struct node

/** getter for p */
inline Pt::pt3D get_p(Node const &n /**< [in] */) { return n.p; }

/** getter for u0*/
inline const Pt::pt3D get_u0(Node const &n /**< [in] */) { return n.u0; }

/** getter for v0*/
inline const Pt::pt3D get_v0(Node const &n /**< [in] */) { return n.v0; }

/** getter for u */
inline const Pt::pt3D get_u(Node const &n /**< [in] */) { return n.u; }

/** getter for v */
inline const Pt::pt3D get_v(Node const &n /**< [in] */) { return n.v; }

/** getter for u component */
inline double get_u_comp(Node const &n /**< [in] */, Pt::index idx /**< [in] */)
    {
    return n.u(idx);
    }

/** getter for v component */
inline double get_v_comp(Node const &n /**< [in] */, Pt::index idx /**< [in] */)
    {
    return n.v(idx);
    }

/** getter for v0 component */
inline double get_v0_comp(Node const &n /**< [in] */, Pt::index idx /**< [in] */)
    {
    return n.v0(idx);
    }

/** getter for phi */
inline double get_phi(Node const &n /**< [in] */) { return n.phi; }

/** getter for phi0 */
inline double get_phi0(Node const &n /**< [in] */) { return n.phi0; }

/** getter for phiv0 */
inline double get_phiv0(Node const &n /**< [in] */) { return n.phiv0; }

/** setter for phi */
inline void set_phi(Node &n, double val) { n.phi = val; }

/** setter for phi_v */
inline void set_phiv(Node &n, double val) { n.phiv = val; }

/** template function to provide P matrix coefficients for T= tetra or facette, with respect to its
 * block diagonal structure */
template<class T>
double Pcoeff(T &x, const int i, const int j)
    {
    const int N = x.getN();
    double val(0);
    int node_i = i % N;

    if (node_i == (j % N))
        {
        const Nodes::Node &n = x.getNode(node_i);

        if (i < N)
            {
            val = n.ep(j / N);
            }
        else
            {
            val = n.eq(j / N);
            }
        }
    return val;
    }

/** template to make projection for T, tetra or facette. It computes Bp = P*B and stores result in
 * inner vector Lp of class T */
template<class T, int N>
void projection_vect(T &x, Pt::pt3D *B)
    {
    for (int i = 0; i < (2 * N); i++)
        {
        x.Lp[i] = 0;
        for (int k = 0; k < N; k++)
            {
            x.Lp[i] += Pcoeff<T>(x,i,k) * B[k].x()
                       + Pcoeff<T>(x,i,N + k) * B[k].y()
                       + Pcoeff<T>(x,i,2*N + k) * B[k].z();
            }
        }
    }

/** template to make projection for T, tetra or facette. It computes Ap = (P*A)*trans(P) and stores
 * result in inner matrix Kp of class T */
template<class T, int N>
void projection_mat(T &x, double (&A)[3 * N][3 * N])
    {
    double PA[2 * N][3 * N];  // no need to initialize with zeros
    for (int i = 0; i < (2 * N); i++)
        {
        for (int k = 0; k < (3 * N); k++)
            {
            PA[i][k] = 0;
            for (int j = 0; j < (3 * N); j++)
                {
                PA[i][k] += Pcoeff<T>(x, i, j) * A[j][k];
                }
            }
        }

    for (int i = 0; i < (2 * N); i++)
        for (int k = 0; k < (2 * N); k++)
            {
            x.Kp[i][k] = 0;
            for (int j = 0; j < (3 * N); j++)
                {
                x.Kp[i][k] += PA[i][j] * Pcoeff<T>(x, k, j);
                }
            }
    }

    }  // end namespace Nodes

#endif /* node_h */
