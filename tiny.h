/** \file tiny.h
contains a namespace to perform simple printing and algebra simple operations on matrices and vectors. <br>
\todo this namespace is redondant with gmm, could probably be replaced or improved with overloaded operators and a matrix class, or valarrays
*/
#include <iostream>
#include <complex>

/**
\namespace tiny to grab altogether the templates to compute some linear algebra
*/
namespace tiny{

/** printing function
*/
template <typename T, int M, int N> inline void print(T A[M][N]) {
   for (int i=0; i<M; i++)
   for (int j=0; j<N; j++) {
       std::cout << A[i][j];
       }
   std::cout << std::endl;
   }

/** scalar product of two vectors
*/
template <typename T, int N> inline T sp(T X[N], T Y[N]) {
   T v=T(0);
   for (int i=0; i<N; i++) 
       v+= X[i]*Y[i];
   return v;
   }
   

/** in place addition += 
\return return in Y \f$ Y \mathrel{+}= X \f$
*/
template <typename T, int N> inline void add(T X[N], T Y[N]) {
   for (int i=0; i<N; i++) 
       Y[i]+= X[i];
   }

/** mat vector multiplication
\return returns Y \f$ Y = A X \f$
*/
template <typename T, int M, int N> inline void mult(T A[M][N], T X[N], T Y[M]) {
   for (int i=0; i<M; i++) {
       T v=T(0);
       for (int j=0; j<N; j++) 
           v+= A[i][j]*X[j];
       Y[i]=v;
       }
   }

/** multiplication of a vector and a transposed matrix
\return returns in Y \f$ Y = A^{\dagger} X \f$
*/
template <typename T, int M, int N> inline void transposed_mult(T X[M], T A[M][N], T Y[N]) {
   for (int j=0; j<N; j++) { 
       T v=T(0);
       for (int i=0; i<M; i++)
           v+= X[i]*A[i][j];
       Y[j]=v;
       }
   }

/** mat mat multiplication <br>
\return \f$ C = A B \f$
*/
template <typename T, int M, int N, int P> inline void mult(T A[M][N], T B[N][P], T C[M][P]) {
   for (int i=0; i<M; i++) 
   for (int k=0; k<P; k++) {
       T v=T(0);
       for (int j=0; j<N; j++) 
           v+= A[i][j]*B[j][k];
       C[i][k]=v;
       }
   }
}

/** 
template to compute average of either u or v on the whole set of nodes
*/
template <int UorV>
double moy(Fem &fem, int d)
{
const int TET = fem.TET;
const int N   = Tet::N;
const int NPI = Tet::NPI;

double sum = 0.;
for (int t=0; t<TET; t++){
    Tet &tet = fem.tet[t];
    double val_nod[N], val[NPI];
    for (int ie=0; ie<N; ie++) {
        int i = tet.ind[ie];
        Node &node = fem.node[i];
	if(UorV)        
		val_nod[ie] = node.u[d];
	else val_nod[ie] = node.v[d]; 
        }
   tiny::transposed_mult<double, N, NPI> (val_nod, tet.a, val);
   sum += tiny::sp<double, NPI> (val, tet.weight);
   }

return sum/fem.vol;
}

