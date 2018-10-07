/** \file linear_algebra.h
\brief secondary header, it grabs altogether the linear algebra by the solver to apply fem method  <br>
It encapsulates the calls to GMM , the assemblage and projection of the matrix for all elements <br>
two templates projection and assemblage template class parameter is either Facette::Fac or Tetra::Tet 
*/

#include "boost/numeric/mtl/mtl.hpp"
#include "boost/numeric/itl/itl.hpp"

#include <boost/numeric/mtl/matrix/inserter.hpp>
#include <boost/numeric/mtl/operation/set_to_zero.hpp>
#include <boost/numeric/mtl/interface/vpt.hpp>

#include "fem.h"

#ifndef linear_algebra_h
#define linear_algebra_h

typedef typename mtl::Collection< mtl::compressed2D<double> >::value_type v_type;
typedef mtl::mat::inserter< mtl::compressed2D<double>,mtl::update_plus<v_type> > sparseInserter;


/** \class LinAlgebra
convenient class to grab altogether some part of the calculations involved using gmm solver at each timestep
*/
class LinAlgebra
{
public:
	/** constructor */	
	inline LinAlgebra(Settings & s,
                      std::vector<Node> & myNode,
                      std::vector <Tetra::Tet> & myTet,
                      std::vector <Facette::Fac> & myFace) 
    {settings = &s; refNode = &myNode; refTet = &myTet; refFac = &myFace;}
	
	//itl::pc::diagonal < mtl::compressed2D<double> > prc;/**< diagonal preconditionner */

	int  vsolve(double dt,long nt);/**< solver */

    /** setter for DW_dz */
    inline void set_DW_vz(double vz){DW_vz = vz;}    

    /** setter for Hext */
    inline void set_Hext(double Hx,double Hy,double Hz){Hext[0]=Hx;Hext[1]=Hy;Hext[2]=Hz;}

    /** getter for v_max */
    inline double get_v_max() {return v_max;}
    
    /** getter node */
    inline Node getNode(int i) {return (*refNode)[i];}
    
    /** getter node physical position */
    inline Pt::pt3D getNodePhysPos(int i) {return (*refNode)[i].p;} 
    
    
private:
    const int MAXITER = 500;/**< maximum number of iteration for biconjugate gradient algorithm */
    const int REFRESH_PRC = 20;/**< refresh every REFRESH_PRC the diagonal preconditioner */
    
    std::vector<Node>  *refNode;/**< direct access to the Nodes */
	std::vector <Facette::Fac> *refFac; /**< direct access to the faces */
	std::vector <Tetra::Tet> *refTet; /**< direct access to the tetrahedrons */
	
	double Hext[DIM];/**< applied field */
    double DW_vz;/**< speed of the domain wall */
	Settings *settings;/**< copy of the settings */
    double v_max;/**< maximum speed */
	
/** computes the local vector basis {ep,eq} in the tangeant plane for projection on the elements */
inline void base_projection(void)
	{ std::for_each(refNode->begin(),refNode->end(),[](Node &n) { n.buildBase_epeq();}); }
	
/**
template function to compute projection of an element <br>
template parameter T is either tetra of face
*/
template <class T>
void projection(T &elt,
           mtl::dense2D <double> &A,  mtl::dense_vector <double> &B,
           mtl::dense2D <double> &Ap, mtl::dense_vector <double> &Bp)
{
const int N = elt.getN();
mtl::dense2D <double> P(2*N,3*N), PA(2*N,3*N);
mtl::mat::set_to_zero(P); mtl::mat::set_to_zero(PA);

for (int i=0; i<N; i++){
    Node &n = (*refNode)[elt.ind[i]];
    P(i,i)  = n.ep.x();  P(i,N+i)  = n.ep.y();  P(i,2*N+i)  = n.ep.z();
    P(N+i,i)= n.eq.x();  P(N+i,N+i)= n.eq.y();  P(N+i,2*N+i)= n.eq.z();
    }

PA = P*A;
Ap = PA*trans(P);
Bp = P*B;
}

/**
template function to perform the
matrix assembly with all the contributions of the tetrahedrons/faces <br>
template parameter T is either tetra or face
*/

template <class T>
void assemblage(T &elt,sparseInserter *ins,
           mtl::dense2D <double> &Ke, mtl::dense_vector <double> &Le,
           mtl::compressed2D<double> &K, mtl::dense_vector<double> &L)
    {
    const int NOD = refNode->size();
    const int N = elt.getN();

    for (int i=0; i < N; i++){
        int i_= elt.ind[i];             
        for (int j=0; j < N; j++){
            int j_= elt.ind[j];
            (*ins)(NOD+i_,j_) << Ke(i,j);      (*ins)(NOD+i_, NOD+j_) << Ke(  i,N+j);
            (*ins)(    i_,j_) << Ke(N+i,j);    (*ins)(    i_, NOD+j_) << Ke(N+i,N+j);
            }
        L(NOD+i_) += Le(i);//L[NOD+i_]+= Le[  i];
        L(i_) += Le(N+i);//L[    i_]+= Le[N+i];
        }
    //std::cout<<"temps assemblage: "<<diff_t<<endl;
    }

    
}; // fin class linAlgebra

#endif
