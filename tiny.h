#ifndef tiny_h
#define tiny_h

/** \file tiny.h
contains a namespace to perform simple printing and algebra simple operations on matrices and vectors. <br>
\todo this namespace is redondant with gmm, could probably be replaced or improved with overloaded operators and a matrix class, or valarrays
*/
#include <iostream>

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

/** mat vector multiplication
\return returns Y \f$ Y = A X \f$
*/
template <typename T, int M, int N> inline void mult(const T A[M][N],const T X[N], T Y[M]) {
   for (int i=0; i<M; i++) {
       T v=T(0);
       for (int j=0; j<N; j++) 
           v+= A[i][j]*X[j];
       Y[i]=v;
       }
   }

/** left multiplication of a transposed vector and a matrix
\return returns in Y \f$ Y = X^{\dagger} A \f$
*/
template <typename T, int M, int N> inline void transposed_mult(const T X[M],const T A[M][N], T Y[N]) {
   for (int j=0; j<N; j++) { 
       T v=T(0);
       for (int i=0; i<M; i++)
           v+= X[i]*A[i][j];
       Y[j]=v;
       }
   }

/** operator: do left multiplication of a transposed vector and a matrix, and multiply by -1
\return returns in Y \f$ Y = -X^{\dagger} A \f$
*/
template <typename T, int M, int N> inline void neg_transposed_mult(const T X[M],const T A[M][N], T Y[N]) {
   for (int j=0; j<N; j++) { 
       T v=T(0);
       for (int i=0; i<M; i++)
           v+= X[i]*A[i][j];
       Y[j] = -v;
       }
   }
   
   
/** const mat const mat multiplication <br>
\return \f$ C = A B \f$
*/
template <typename T, int M, int N, int P> inline void mult(const T A[M][N],const T B[N][P], T C[M][P]) {
   for (int i=0; i<M; i++) 
   for (int k=0; k<P; k++) {
       T v=T(0);
       for (int j=0; j<N; j++) 
           v+= A[i][j]*B[j][k];
       C[i][k]=v;
       }
   }
}

#endif
