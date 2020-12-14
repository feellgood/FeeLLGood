#ifndef tiny_h
#define tiny_h

/** \file tiny.h
contains a namespace to perform simple printing and algebra simple operations on matrices and vectors. <br>
\todo this namespace is redondant with gmm, could probably be replaced or improved with overloaded operators and a matrix class, or valarrays

addition, substraction, frobenius norms and dist implemented only for unit tests, not used in feellgood algebra
*/
#include <iostream>

#include "pt3D.h"

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
   
/** frobenius Norm for vector with real coefficients
*/
template <typename T,int N> inline double frob_norm(const T A[N])
{
T result=T(0);

for (int j=0; j<N; j++) 
    result += Pt::sq(A[j]);

return sqrt(result);
}   

/** frobenius Norm for matrix with real coefficients
*/
template <typename T,int M,int N> inline double frob_norm(const T A[M][N])
{
T result=T(0);

for (int i=0; i<M; i++)
    {
    for (int j=0; j<N; j++) 
        result += Pt::sq(A[i][j]);
       }
return sqrt(result);
}

/** distance defined with frobenius Norm for matrix with real coefficients : dist(A,B) = Frob(A-B)
*/
template <typename T,int M,int N> inline double dist(const T A[M][N],const T B[M][N])
{
T result=T(0);

for (int i=0; i<M; i++)
    for (int j=0; j<N; j++) 
        result += Pt::sq(A[i][j] - B[i][j]);

return sqrt(result);
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
   
/** mat Transpose(mat) multiplication
\return \f$ C = A B^{\dagger}  \f$
*/
template <typename T, int M, int N, int P> inline void direct_transposed_mult(T A[M][N],T B[P][N], T C[M][P]) {
   for (int i=0; i<M; i++) 
   for (int k=0; k<P; k++) {
       T v=T(0);
       for (int j=0; j<N; j++) 
           v+= A[i][j]*B[k][j];
       C[i][k]=v;
       }
   }

   /** mat mat multiplication
\return \f$ C = A B \f$
*/
template <typename T, int M, int N, int P> inline void mult(const Pt::pt3D vec[N],const T B[N][P], Pt::pt3D C[P]) {
   for (int i=0; i<M; i++) 
   for (int k=0; k<P; k++) {
       T v=T(0);
       for (int j=0; j<N; j++) 
           v+= vec[j](i)*B[j][k];
       C[k](i,v);
       }
   }
   
/** mat mat multiplication
\return \f$ C = A B \f$
*/
template <typename T, int M, int N, int P> inline void mult(const Pt::pt3D vec[N],const T B[N][P], T C[M][P]) {
   for (int i=0; i<M; i++) 
   for (int k=0; k<P; k++) {
       T v=T(0);
       for (int j=0; j<N; j++) 
           v+= vec[j](i)*B[j][k];
       C[i][k]=v;
       }
   }

/** mat mat addition
\return \f$ C = A + B \f$
*/
template <typename T, int M, int N> inline void add(const T A[M][N],const T B[M][N], T C[M][N])
{
for (int i=0; i<M; i++) 
   for (int k=0; k<N; k++) 
    { C[i][k] = A[i][k] + B[i][k]; }
}

/** mat mat substraction
\return \f$ C = A - B \f$
*/
template <typename T, int M, int N> inline void sub(const T A[M][N],const T B[M][N], T C[M][N])
{
for (int i=0; i<M; i++) 
   for (int k=0; k<N; k++) 
    { C[i][k] = A[i][k] - B[i][k]; }
}

/** mat mat multiplication
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
