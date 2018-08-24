/** \file linear_algebra.h
\brief secondary header, it grabs altogether the linear algebra by the solver to apply fem method  <br>
It encapsulates the calls to GMM , the assemblage and projection of the matrix for all elements <br>
two templates projection and assemblage template class parameter is either Facette::Fac or Tetra::Tet 
*/

#include "gmm/gmm_kernel.h"  // ct gmm_kernel.h plutôt que gmm.h , qui appelle des fichiers de getfem
#include "gmm/gmm_precond_diagonal.h" //ct

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
	inline LinAlgebra(Fem &f, Settings &s) {fem = f; settings = s;}
	gmm::diagonal_precond <read_matrix>  *prc;/**< diagonal preconditionner */

int  vsolve(long nt);/**< solver */
void base_projection(void);/**< computes the vector basis projection on the elements */

private:
	Fem fem;/**< access to some part of struct fem */
	Settings settings;/**< copy of the settings */
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
    Node &n = fem.node[elt.ind[i]];
    P(i,i)  = n.ep[0];  P(i,N+i)  = n.ep[1];  P(i,2*N+i)  = n.ep[2];
    P(N+i,i)= n.eq[0];  P(N+i,N+i)= n.eq[1];  P(N+i,2*N+i)= n.eq[2];
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
const int NOD = fem.NOD;
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
