/** \file linear_algebra.h
\brief secondary header, it grabs altogether the linear algebra by the solver to apply fem method  <br>
It encapsulates the calls to GMM , the assemblage and projection of the matrix for all elements <br>
two templates projection and assemblage template class parameter is either Facette::Fac or Tetra::Tet 
*/

#include "gmm_kernel.h"  // ct gmm_kernel.h plutôt que gmm.h , qui appelle des fichiers de getfem
#include "gmm_precond_diagonal.h" //ct

#include "fem.h"

#ifndef linear_algebra_h
#define linear_algebra_h

typedef gmm::wsvector <double>   write_vector;/**< convenient typedef : a write sparse vector in GMM is a stl::map, complexity for both read and write an element is \f$ \log(N) \f$ */
typedef gmm::rsvector <double>   read_vector;/**< convenient typedef : a read sparse vector */

typedef gmm::row_matrix	<write_vector>   write_matrix;/**< convenient typedef for row matrix in write mode */
typedef gmm::row_matrix	<read_vector>    read_matrix;/**< convenient typedef for row matrix in read mode */

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
	
	gmm::diagonal_precond <read_matrix>  * prc;/**< diagonal preconditionner */

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
           gmm::dense_matrix <double> &A,  std::vector <double> &B,
           gmm::dense_matrix <double> &Ap, std::vector <double> &Bp)
{
const int N = elt.getN();
gmm::dense_matrix <double> P(2*N,3*N), PA(2*N,3*N);
for (int i=0; i<N; i++){
    Node &n = (*refNode)[elt.ind[i]];
    P(i,i)  = n.ep.x();  P(i,N+i)  = n.ep.y();  P(i,2*N+i)  = n.ep.z();
    P(N+i,i)= n.eq.x();  P(N+i,N+i)= n.eq.y();  P(N+i,2*N+i)= n.eq.z();
    }

mult(P,A,PA);
mult(PA, gmm::transposed(P), Ap);

mult(P,B,Bp);
}

/**
template function to perform the
matrix assembly with all the contributions of the tetrahedrons/faces <br>
template parameter T is either tetra or face
*/
template <class T>
void assemblage(T &elt,
           gmm::dense_matrix <double> &Ke, std::vector <double> &Le,
           write_matrix &K, write_vector &L)
    {
    const int NOD = refNode->size();
    const int N = elt.getN();

    for (int i=0; i < N; i++){
        int i_= elt.ind[i];             
        for (int j=0; j < N; j++){
            int j_= elt.ind[j];
            K(NOD+i_,j_)+= Ke(i,j);      K(NOD+i_, NOD+j_)+= Ke(  i,N+j);
            K(    i_,j_)+= Ke(N+i,j);    K(    i_, NOD+j_)+= Ke(N+i,N+j);
            }
        L[NOD+i_]+= Le[  i];
        L[    i_]+= Le[N+i];
        }
    //std::cout<<"temps assemblage: "<<diff_t<<endl;
    }

    
}; // fin class linAlgebra

#endif
