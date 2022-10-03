#ifndef node_h
#define node_h

/** \file node.h
\brief header to define struct Node 
*/

#include <memory>

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

Pt::pt3D ep;/**< local vector basis : \f$ e_p = \vec{rand} \times u0 \f$ , then normalized */
Pt::pt3D eq;/**< local vector basis : \f$ e_q = u0 \times e_p \f$ , then normalized */

double phi0;/**< scalar potential initial or reset value, used to store previous value for time evolution */
double phi;/**< scalar potential value */
double phiv0;/**< initial or reset value, used to store previous value for time evolution */
double phiv;/**< scalar potential associated to velocity */

double V;/**< electrostatic potential (for STT) */

/**
computes vector ep, it is the second vector of a base composed of  \f$ \left\lbrace u0,e_p , e_q = u0 \times e_p \right\rbrace \f$
 */
inline void calc_ep() {ep = Pt::pt3D(theta_sph,phi_sph)*u0; ep.normalize(); }

/**
computes vector eq, it is the third vector of a base composed of  \f$ \left\lbrace u0,e_p , e_q = u0 \times e_p \right\rbrace \f$
 */
inline void calc_eq() {eq = u0*ep; ep.normalize(); }

/** setter for the local basis vector */
inline void setBasis(const double r1,const double r2)
{
theta_sph = r1;
phi_sph = r2;
calc_ep();
calc_eq();
}

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
in a base composed of u0,ep,eq = u0*ep we have  
\f$ v = v_p e_p + v_q e_q \f$
and new magnetization value is : \f$ u = u_0 + v dt \f$ after normalization
There is an algebraic simplification when developping \f$ v \f$ due to the fact that magnetization u0 is a unit vector, and \f$ v = v_p/r \hat{\zeta} \times u_0 + v_q/r (\hat{\zeta} - (u_0 \cdot \hat{\zeta}) u_0 )  \f$
With \f$ \hat{\zeta} \f$ unit vector defined by theta_sph and phi_sph and \f$ r=sin(\hat{\zeta},u_0) \f$
This formula needs less operator*
*/

inline void make_evol(const double vp /**< [in] */,const double vq /**< [in] */,const double dt /**< [in] */)
{
Pt::pt3D zeta = Pt::pt3D(theta_sph,phi_sph); 
double r = Pt::norme(zeta*u0);

v = vp*ep + (vq/r)*(zeta - (Pt::pScal(u0,zeta))*u0);
u = u0 + dt*v;
u.normalize();
}

};//end struct node

/** getter for p */
inline Pt::pt3D get_p(Node const& n /**< [in] */) {return n.p;}

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

/** template function to provide P matrix coefficients for T= tetra or facette, with respect to its block diagonal structure */
template<class T> double Pcoeff(T & x,const int i,const int j)
    {
    const int N = x.getN();
    double val(0);
    int node_i = i%N;
    
    if(node_i == (j%N))
        {
        const Nodes::Node & n = x.getNode(node_i);
           
        if(i<N)
            { val = n.ep(j/N); }
        else
            { val = n.eq(j/N); }
        }
    return val;
    }

/** template to make projection for T, tetra or facette. It computes Bp = P*B and stores result in inner vector Lp of class T*/
    //template<class T,int N> void projection_vect(T &x, Pt::pt3D (&B)[N])
template<class T,int N> void projection_vect(T &x, Pt::pt3D *B)
    {
        //tiny::mult<double,2*N,3*N>(P,B,x.Lp);
    for (int i=0; i<(2*N); i++)
        {
        x.Lp[i] = 0;
        for (int k=0; k<N; k++) { x.Lp[i] += Nodes::Pcoeff<T>(x,i,k)*B[k].x(); }
        for (int k=0; k<N; k++) { x.Lp[i] += Nodes::Pcoeff<T>(x,i,N+k)*B[k].y(); }
        for (int k=0; k<N; k++) { x.Lp[i] += Nodes::Pcoeff<T>(x,i,2*N+k)*B[k].z(); }
        }
    }
    
    
/** template to make projection for T, tetra or facette. It computes Ap = (P*A)*trans(P) and stores result in inner matrix Kp of class T*/
template <class T,int N> void projection_mat(T &x,double (&A)[3*N][3*N])
	{
    double PA[2*N][3*N]; // no need to initialize with zeros
//tiny::mult<double,2*N,3*N,3*N>(P,A,PA);
	for (int i=0; i<(2*N); i++) 
        {
        for (int k=0; k<(3*N); k++)
            {
            PA[i][k]=0;
            for (int j=0; j<(3*N); j++) { PA[i][k] += Nodes::Pcoeff<T>(x,i,j)*A[j][k]; }
            }
        }

//tiny::direct_transposed_mult<double,2*N,3*N,2*N>(PA,P,x.Kp);
   for (int i=0; i<(2*N); i++) 
        for (int k=0; k<(2*N); k++)
        {
       x.Kp[i][k] = 0;
       for (int j=0; j<(3*N); j++) { x.Kp[i][k] += PA[i][j]* Nodes::Pcoeff<T>(x,k,j); }
       }
    }

/**
 * This is a std::shared_ptr<Node> with the addition of operator[],
 * for compatibility with c++14 (operator[] is defined in C++17).
 */
class NodeList : public std::shared_ptr<Node>
{
public:
    /**
     * Constructor from a regular pointer.
     *
     * On C++17, we could use std::shared_ptr<Node[]>, which would use
     * the appropriate array deleter (delete[]). Since we want
     * compatibility with C++14, we have to use std::shared_ptr<Node>,
     * and specify a custom deleter, otherwise the "scalar" delete would
     * be used.
     */
    NodeList(Node *p = nullptr)
    : std::shared_ptr<Node>(p, std::default_delete<Node[]>())
    {}

    /** Array subscript operator */
    Node& operator[](std::ptrdiff_t i) { return (&**this)[i]; }

    /** Array subscript operator, `const` version */
    const Node& operator[](std::ptrdiff_t i) const { return (&**this)[i]; }
};
    
} //end namespace Nodes

#endif /* node_h */
