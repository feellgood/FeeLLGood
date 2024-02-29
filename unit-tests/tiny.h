#ifndef tiny_h
#define tiny_h

/** \file tiny.h
contains a namespace to perform simple printing and algebra simple operations on matrices and
vectors. <br> \todo this namespace is redondant with eigen, could probably be replaced or improved
with overloaded operators and a matrix class, or valarrays

addition, substraction, frobenius norms and dist implemented only for unit tests, not used in
feellgood algebra
*/
#include <iostream>
#include <cmath>

double _sq(const double x) { return x * x; }

/**
\namespace tiny to grab altogether the templates to compute some linear algebra
*/
namespace tiny
    {

/** printing function
 */
template<typename T, int M, int N>
inline void print(T A[M][N])
    {
    for (int i = 0; i < M; i++)
        for (int j = 0; j < N; j++)
            {
            std::cout << A[i][j];
            }
    std::cout << std::endl;
    }

/** frobenius Norm for vector with real coefficients
 */
template<typename T, int N>
inline double frob_norm(const T A[N])
    {
    T result = T(0);

    for (int j = 0; j < N; j++)
        result += _sq(A[j]);

    return sqrt(result);
    }

/** frobenius Norm for matrix with real coefficients
 */
template<typename T, int M, int N>
inline double frob_norm(const T A[M][N])
    {
    T result = T(0);

    for (int i = 0; i < M; i++)
        {
        for (int j = 0; j < N; j++)
            result += _sq(A[i][j]);
        }
    return sqrt(result);
    }

/** distance defined with frobenius Norm for matrix with real coefficients : dist(A,B) = Frob(A-B)
 */
template<typename T, int M, int N>
inline double dist(const T A[M][N], const T B[M][N])
    {
    T result = T(0);

    for (int i = 0; i < M; i++)
        for (int j = 0; j < N; j++)
            result += _sq(A[i][j] - B[i][j]);

    return sqrt(result);
    }

/** matrix-vector multiplication: \f$ Y = A X \f$
\param[in] A input matrix
\param[in] X input vector
\param[out] Y result
*/
template<typename T, int M, int N>
inline void mult(const T A[M][N], const T X[N], T Y[M])
    {
    for (int i = 0; i < M; i++)
        {
        T v = T(0);
        for (int j = 0; j < N; j++)
            v += A[i][j] * X[j];
        Y[i] = v;
        }
    }

/** left multiplication of a transposed vector and a matrix: \f$ Y = X^{\dagger} A \f$
\param[in] X input vector
\param[in] A input matrix
\param[out] Y result
*/
template<typename T, int M, int N>
inline void transposed_mult(const T X[M], const T A[M][N], T Y[N])
    {
    for (int j = 0; j < N; j++)
        {
        T v = T(0);
        for (int i = 0; i < M; i++)
            v += X[i] * A[i][j];
        Y[j] = v;
        }
    }

/** left multiplication of a transposed vector and a matrix, multiplied by -1:
    \f$ Y = -X^{\dagger} A \f$
\param[in] X input vector
\param[in] A input matrix
\param[out] Y result
*/
template<typename T, int M, int N>
inline void neg_transposed_mult(const T X[M], const T A[M][N], T Y[N])
    {
    for (int j = 0; j < N; j++)
        {
        T v = T(0);
        for (int i = 0; i < M; i++)
            v += X[i] * A[i][j];
        Y[j] = -v;
        }
    }

/** multiplication of a matrix with the transpose of another matrix:
    \f$ C = A B^{\dagger}  \f$
\param[in] A,B input matrices
\param[out] C result
*/
template<typename T, int M, int N, int P>
inline void direct_transposed_mult(T A[M][N], T B[P][N], T C[M][P])
    {
    for (int i = 0; i < M; i++)
        for (int k = 0; k < P; k++)
            {
            T v = T(0);
            for (int j = 0; j < N; j++)
                v += A[i][j] * B[k][j];
            C[i][k] = v;
            }
    }

/** matrix addition: \f$ C = A + B \f$
\param[in] A,B input matrices
\param[out] C result
*/
template<typename T, int M, int N>
inline void add(const T A[M][N], const T B[M][N], T C[M][N])
    {
    for (int i = 0; i < M; i++)
        for (int k = 0; k < N; k++)
            {
            C[i][k] = A[i][k] + B[i][k];
            }
    }

/** matrix subtraction: \f$ C = A - B \f$
\param[in] A,B input matrices
\param[out] C result
*/
template<typename T, int M, int N>
inline void sub(const T A[M][N], const T B[M][N], T C[M][N])
    {
    for (int i = 0; i < M; i++)
        for (int k = 0; k < N; k++)
            {
            C[i][k] = A[i][k] - B[i][k];
            }
    }

/** matrix multiplication: \f$ C = A B \f$
\param[in] A,B input matrices
\param[out] C result
*/
template<typename T, int M, int N, int P>
inline void mult(const T A[M][N], const T B[N][P], T C[M][P])
    {
    for (int i = 0; i < M; i++)
        for (int k = 0; k < P; k++)
            {
            T v = T(0);
            for (int j = 0; j < N; j++)
                v += A[i][j] * B[j][k];
            C[i][k] = v;
            }
    }
    }  // namespace tiny

#endif
