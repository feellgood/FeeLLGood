#ifndef node_h
#define node_h

/** \file node.h
\brief header to define struct Node 
*/

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
        {IDX_p = 0,IDX_u0 = 1,IDX_v0 = 2,IDX_u = 3,IDX_v = 4,IDX_theta_sph = 5,IDX_phi_sph = 6,IDX_phi0 = 7,IDX_phi = 8,IDX_phiv0 = 9,IDX_phiv = 10};

/** \struct Node
Node is containing physical point of coordinates \f$ p = (x,y,z) \f$, magnetization value at \f$ m(p,t) \f$. 
Many other values for the computation of the scalar potential \f$ \phi \f$
Angles theta_sph and phi_sph are defining a unit vector to build a local base u0,ep = u0*e(theta_sph,phi_sph), eq = u0*ep
*/
struct Node {
Pt::pt3D p;/**< Physical position p=(x,y,z)  of the node */
Pt::pt3D u0;/**< magnetization initial or reset value, used to store previous value for time evolution */
Pt::pt3D v0;/**< initial or reset value, used to store previous value for time evolution */
Pt::pt3D u;/**< magnetization value */
Pt::pt3D v;/**< magnetization speed */

double theta_sph;/**< theta angle for unit vector \f$ e_p = (\theta,\phi) \times u_0 \f$ in spherical coordinates  */
double phi_sph;/**< phi angle for unit vector \f$ e_p = (\theta,\phi) \times u_0 \f$ in spherical coordinates  */

double phi0;/**< scalar potential initial or reset value, used to store previous value for time evolution */
double phi;/**< scalar potential value */
double phiv0;/**< initial or reset value, used to store previous value for time evolution */
double phiv;/**< no idea */

/**
vector ep is computed when needed, it is the second vector of a base composed of u0,ep,eq = u0*ep , vector product 
 */
inline Pt::pt3D calc_ep() const {Pt::pt3D ep = Pt::pt3D(theta_sph,phi_sph)*u0; ep.normalize(); return ep; }


/**
vector eq is computed when needed, it is the third vector of a base composed of u0,ep,eq = u0*ep , vector product 
 */
inline Pt::pt3D calc_eq() const {return u0*calc_ep(); }

/**
reset the node magnetization, speed, phi, and phiv
*/
inline void reset(void) { u=u0; v=v0; phi = phi0; phiv = phiv0;}

/** 
preparation of the quantities u0,v0,phi0,phiv0 for incomming time-step
*/
inline void evolution(void) { u0=u; v0=v; phi0 = phi; phiv0 = phiv; }

/**
integration of the evolution of the magnetization for time step dt
*/
inline void make_evol(double vp /**< [in] */,double vq /**< [in] */,double dt /**< [in] */) {
	Pt::pt3D ep = calc_ep(); 
    v = vp*ep + vq*u0*ep;
	u = u0 + dt*v;
	u.normalize();}

};//end struct node

/** getter for p */
inline const Pt::pt3D get_p(Node const& n /**< [in] */) {return n.p;}

/** getter for u0*/
inline const Pt::pt3D get_u0(Node const& n /**< [in] */) {return n.u0;}

/** getter for v0*/
inline const Pt::pt3D get_v0(Node const& n /**< [in] */) {return n.v0;}

/** getter for u */
inline const Pt::pt3D get_u(Node const& n /**< [in] */) {return n.u;}

/** getter for v */
inline const Pt::pt3D get_v(Node const& n /**< [in] */) {return n.v;}

/** getter for u component */
inline double get_u_comp(Node const& n /**< [in] */, Pt::index idx /**< [in] */) {return n.u(idx);}

/** getter for v component */
inline double get_v_comp(Node const& n /**< [in] */ ,Pt::index idx /**< [in] */) {return n.v(idx);}


/** getter for phi */
inline double get_phi(Node const& n /**< [in] */) {return n.phi;}

/** getter for phi0 */
inline double get_phi0(Node const& n /**< [in] */) {return n.phi0;}

/** getter for phiv0 */
inline double get_phiv0(Node const& n /**< [in] */) {return n.phiv0;}

/** setter for phi */
inline void set_phi(Node & n,double val) {n.phi = val;}

/** setter for phi_v */
inline void set_phiv(Node & n,double val) {n.phiv = val;}



} //end namespace Nodes

#endif /* node_h */
