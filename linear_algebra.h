#include "gmm/gmm_kernel.h"  // ct gmm_kernel.h plutôt que gmm.h , qui appelle des fichiers de getfem
#include "gmm/gmm_precond_diagonal.h" //ct

#include "fem.h"

#ifndef linear_algebra_h
#define linear_algebra_h

typedef gmm::wsvector <double>   write_vector;/**< convenient macro */
typedef gmm::rsvector <double>   read_vector;/**< convenient macro */

typedef gmm::row_matrix	<write_vector>   write_matrix;/**< convenient macro */
typedef gmm::row_matrix	<read_vector>    read_matrix;/**< convenient macro */

/** \class linAlgebra
convenient class to grab altogether some part of the calculations involved using gmm for the solver at each timestep
*/
class LinAlgebra
{
public:
	inline LinAlgebra() {}
	gmm::diagonal_precond <read_matrix>  *prc;/**< diagonal preconditionner */
	//gmm::ilutp_precond < read_matrix > *prc;

/**
computes the contribution of the tetrahedron to the integrals
*/
void integrales(Fem &fem,Settings &settings, Tet &tet, gmm::dense_matrix <double> &AE, vector <double> &BE);

/**
computes the contribution of the surface to the integrals
*/
void integrales(Fem &fem,Settings &settings, Fac &fac, gmm::dense_matrix <double> &AE, vector <double> &BE);

int  vsolve(Fem &fem,Settings &settings, long nt);/**< solver */
void base_projection(Fem &fem);/**< computes the projection of the llg operators on the elements */

private:
/**
template function to compute projection of an element <br>
template parameter T is either tetra of face
*/
template <class T>
void projection(Fem &fem, T &elt,
           gmm::dense_matrix <double> &A,  vector <double> &B,
           gmm::dense_matrix <double> &Ap, vector <double> &Bp)
{
const int N = T::N;
gmm::dense_matrix <double> P(2*N,3*N), PA(2*N,3*N);
for (int i=0; i<N; i++){
    Node &node = fem.node[elt.ind[i]];
    P(i,i)  = node.ep[0];  P(i,N+i)  = node.ep[1];  P(i,2*N+i)  = node.ep[2];
    P(N+i,i)= node.eq[0];  P(N+i,N+i)= node.eq[1];  P(N+i,2*N+i)= node.eq[2];
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
void assemblage(Fem &fem, T &elt,
           gmm::dense_matrix <double> &Ke, vector <double> &Le,
           write_matrix &K, write_vector &L)
{
const int NOD = fem.NOD;
const int N = T::N;

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
