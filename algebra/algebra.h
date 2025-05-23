#ifndef ALGEBRA_H
#define ALGEBRA_H

/** \file algebra.h
 * \brief set of class to handle sparse matrix operations for gradient conjugate algorithm
 * a sparse vector class
 * a read and a write sparse matrix class
 * various vector operations, scalar and direct products; ...
 * most of the member functions of the classes are using lambdas and C++11/17 algorithm and numeric
 * */

#include <iostream>
#include <cmath> // sqrt,fabs
#include <algorithm>


#include "coeff.h"
#include "sparseVect.h"
#include "sparseMat.h"

#include "algebraCore.h"
#include "iter.h"

/** \namespace algebra
 * grab altogether sparse matrix, vector and dedicated functions, and cg and bicgstab algorithms
*/

namespace algebra
{
/** returns x² */
template <typename T>
T sq(const T x) { return x * x; }

/** Y *= alpha */
template <typename T>
void scaled( const T alpha, Vector<T> & Y)
    { std::for_each(Y.begin(),Y.end(),[alpha](T &_x){ _x *= alpha; }); }

/** direct product : Z = X⊗Y */
template <typename T>
void p_direct(const Vector<T> & X,const Vector<T> & Y, Vector<T> & Z)
    { std::transform(X.begin(),X.end(),Y.begin(),Z.begin(), std::multiplies<T>() ); }

/** Y += X */
template <typename T>
void add(const Vector<T> & X, Vector<T> & Y)
    {
    for (size_t i = 0; i < Y.size(); ++i)
        Y[i] += X[i];
    }

/** Y -= X */
template <typename T>
void sub(const Vector<T> & X, Vector<T> & Y)
    {
    for (size_t i = 0; i < Y.size(); ++i)
        Y[i] -= X[i];
    }

/** Y += alpha*X       */
template <typename T>
void scaled_add(const Vector<T> & X,const T alpha, Vector<T> & Y)
    {
    for (size_t i = 0; i < Y.size(); ++i)
        Y[i] += alpha * X[i];
    }

/** Y = A*X with r_sparseMat A */
template <typename T>
void mult(r_sparseMat & A, Vector<T> const& X, Vector<T> &Y)
    {
    const int _size = X.size();
    Y.resize(_size);
    if (A.getDim() == _size)
        { for(int i=0;i<_size;i++) Y[i]= A(i).dot(X); }
    }

/** apply a mask to vector X : all coefficients in vector mask are zeroed in X */
template <typename T>
void applyMask(const Vector<int>& mask, Vector<T> & X)
    { std::for_each(mask.begin(),mask.end(),[&X](const int _i){ X[_i] = (T)(0); }); }

/** generic template for gradient conjugate */
template<bool MASK,typename T>
T generic_cg(iteration<T> &iter, r_sparseMat& A, Vector<T> & x, const Vector<T> & rhs, const Vector<T> & xd, const std::vector<int>& ld);

/** conjugate gradient with diagonal preconditioner, returns residu. Template generic_cg is called with dummy xd and ld */
template <typename T>
T cg(iteration<T> &iter, r_sparseMat& A, Vector<T> & x, Vector<T> & b)
    {
    Vector<T> xd;
    std::vector<int> ld;
    return generic_cg<false,T>(iter,A,x,b,xd,ld);
    }

/** conjugate gradient with diagonal preconditioner with Dirichlet conditions, returns residu */
template <typename T>
T cg_dir(iteration<T> &iter, r_sparseMat& A, Vector<T> & x, const Vector<T> & rhs,
         const Vector<T> & xd, const std::vector<int>& ld)
    { return generic_cg<true,T>(iter,A,x,rhs,xd,ld); }

/** operator<< for r_sparseVect */
inline std::ostream & operator<<(std::ostream & flux, r_sparseVect const& v)
    {v.print(flux); return flux;}

/** operator<< for r_sparseMat */
inline std::ostream & operator<<(std::ostream & flux, r_sparseMat const& m)
    {m.print(flux); return flux;}

template <typename T>
inline bool check(Vector<T> &v)
    { return std::none_of(v.begin(),v.end(), [](T &x){ return std::isnan(x);} );  }

} // end namespace algebra

#endif

