/** \file node.h
\brief header to define struct Node 
*/

#ifndef node_h
#define node_h

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
        {IDX_p = 0,IDX_u0 = 1,IDX_v0 = 2,IDX_u = 3,IDX_v = 4,IDX_ep = 5,IDX_eq = 6,IDX_phi0 = 7,IDX_phi = 8,IDX_phiv0 = 9,IDX_phiv = 10};

/** \struct Node
Node is containing physical point of coordinates \f$ p = (x,y,z) \f$, magnetization value at \f$ m(p,t) \f$. 
Many other values for the computation of the scalar potential \f$ \phi \f$
*/
struct Node {
Pt::pt3D p;/**< Physical position p=(x,y,z)  of the node */
Pt::pt3D u0;/**< magnetization initial or reset value, used to store previous value for time evolution */
Pt::pt3D v0;/**< initial or reset value, used to store previous value for time evolution */
Pt::pt3D u;/**< magnetization value */
Pt::pt3D v;/**< magnetization speed */
Pt::pt3D ep;/**< base vector */
Pt::pt3D eq;/**< second base vector */
double phi0;/**< scalar potential initial or reset value, used to store previous value for time evolution */
double phi;/**< scalar potential value */
double phiv0;/**< initial or reset value, used to store previous value for time evolution */
double phiv;/**< no idea */

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
inline void make_evol(double vp,double vq,double dt) {
	v = vp*ep + vq*eq;
	u = u0 + dt*v;
	u.normalize();}

/**
local vector base {ep,eq} in the magnetization tangent plane is built here
*/
inline void buildBase_epeq()
{
if ( Pt::norme2( u0 ) > 0)
	{         
    ep = Pt::rand()*u0;
	ep.normalize();        
	eq = u0*ep;        
	}
    else{ ep = Pt::pt3D(0.,0.,0.); eq = Pt::pt3D(0.,0.,0.);} 
}


};//end struct node

/** getter for u0*/
inline const Pt::pt3D get_u0(Node const& n) {return n.u0;}

/** getter for v0*/
inline const Pt::pt3D get_v0(Node const& n) {return n.v0;}

/** getter for u */
inline const Pt::pt3D get_u(Node const& n) {return n.u;}

/** getter for u component */
inline double get_u_comp(Node n,Pt::index idx) {return n.u(idx);}

/** getter for v component */
inline double get_v_comp(Node n,Pt::index idx) {return n.v(idx);}


/** getter for phi */
inline double get_phi(Node const& n) {return n.phi;}

/** getter for phi0 */
inline double get_phi0(Node const& n) {return n.phi0;}

/** getter for phiv0 */
inline double get_phiv0(Node const& n) {return n.phiv0;}


} //end namespace Nodes

#endif /* node_h */
