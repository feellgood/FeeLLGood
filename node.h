/** \file node.h
\brief header to define struct Node 
*/

#ifndef node_h
#define node_h

#include "pt3D.h"

#include "feellgoodSettings.h"

/** \struct Node
Node is containing physical point of coordinates \f$ p = (x,y,z) \f$, magnetization value at \f$ m(p,t) \f$. 
Many other values for the computation of the scalar potential \f$ \phi \f$
*/
struct Node {
Pt::pt3D p;/**< Physical position p=(x,y,z)  of the node */

triple u0;/**< magnetization initial or reset value, used to store previous value for time evolution */
triple v0;/**< initial or reset value, used to store previous value for time evolution */
triple u;/**< magnetization value */
triple v;/**< magnetization speed */
triple ep;/**< base vector */
triple eq;/**< second base vector */
double phi0;/**< scalar potential initial or reset value, used to store previous value for time evolution */
double phi;/**< scalar potential value */
double phiv0;/**< initial or reset value, used to store previous value for time evolution */
double phiv;/**< no idea */

/**
reset the node magnetization, speed, phi, and phiv
*/
inline void reset(void) {
	u[0] = u0[0]; u[1] = u0[1]; u[2] = u0[2]; 
	v[0] = v0[0]; v[1] = v0[1]; v[2] = v0[2]; 
	phi = phi0; phiv = phiv0;}
};

#endif /* node_h */
