#include <iostream>
#include <complex>

using namespace std;

namespace tiny{

template <typename T, int M, int N> inline void print(T A[M][N]) {
   for (int i=0; i<M; i++)
   for (int j=0; j<N; j++) {
       cout << A[i][j];
       }
   cout << endl;
   }

template <typename T, int N> inline double sp(T X[N], T Y[N]) {
   T v=T(0);
   for (int i=0; i<N; i++) 
       v+= X[i]*Y[i];
   return v;    // FIXME: the return type should be T, right?
   }
   

template <typename T, int N> inline void add(T X[N], T Y[N]) {
   for (int i=0; i<N; i++) 
       Y[i]+= X[i];
   }

template <typename T, int M, int N> inline void mult(T A[M][N], T X[N], T Y[M]) {
   for (int i=0; i<M; i++) {
       T v=T(0);
       for (int j=0; j<N; j++) 
           v+= A[i][j]*X[j];
       Y[i]=v;
       }
   }

template <typename T, int M, int N> inline void transposed_mult(T X[M], T A[M][N], T Y[N]) {
   for (int j=0; j<N; j++) { 
       T v=T(0);
       for (int i=0; i<M; i++)
           v+= X[i]*A[i][j];
       Y[j]=v;
       }
   }

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
