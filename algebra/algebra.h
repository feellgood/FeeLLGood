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
#include <numeric> // inner_product

#include "coeff.h"
#include "sparseVect.h"
#include "sparseMat.h"
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
void scaled( const T alpha, std::vector<T> & Y)
    { std::for_each(Y.begin(),Y.end(),[alpha](T &_x){ _x *= alpha; }); }

/** returns scalar product X.Y */
template <typename T>
auto dot(const std::vector<T,std::allocator<T>> & X,const std::vector<T,std::allocator<T>> & Y)
    {
    auto val = std::inner_product(X.begin(),X.end(),Y.begin(),(T)(0));
    if(!std::isfinite(val)) {std::cout <<"dot:Warning: NaN or inf value.\n"; exit(1);}
    return val;
    }

/** direct product : Z = X⊗Y */
template <typename T>
void p_direct(const std::vector<T> & X,const std::vector<T> & Y,std::vector<T> & Z)
    { for(unsigned int i=0;i<Z.size();i++) Z[i]=X[i]*Y[i]; }

/** Y += X */
template <typename T>
void add(const std::vector<T> & X, std::vector<T> & Y)
    { std::transform(Y.begin(),Y.end(),X.begin(),Y.begin(),std::plus<T>()  ); }

/** Y -= X */
template <typename T>
void sub(const std::vector<T> & X, std::vector<T> & Y)
    { std::transform(Y.begin(),Y.end(),X.begin(),Y.begin(),std::minus<T>()  ); }

/** Y += alpha*X       */
template <typename T>
void scaled_add(const std::vector<T> & X,const T alpha, std::vector<T> & Y)
    { std::transform(Y.begin(),Y.end(),X.begin(),Y.begin(),[alpha] (const T _x,T _y) { return _x+(alpha*_y); }   ); }

/** euclidian norm of vector X */
template <typename T>
auto norm(std::vector<T,std::allocator<T>> & X) { return sqrt(fabs( dot<T>(X,X) )); }

/** Y = A*X with r_sparseMat A */
template <typename T>
void mult(r_sparseMat & A,std::vector<T> const& X,std::vector<T> &Y)
    {
    const int _size = X.size();
    Y.resize(_size);
    if (A.getDim() == _size)
        { for(int i=0;i<_size;i++) Y[i]= A(i).dot(X); }
    }

/** apply a mask to vector X : all coefficients in vector mask are zeroed in X */
template <typename T>
void applyMask(const std::vector<int>& mask, std::vector<T> & X)
    { std::for_each(mask.begin(),mask.end(),[&X](const int _i){ X[_i] = (T)(0); }); }

/** generic template for gradient conjugate */
template<bool MASK,typename T>
T generic_cg(iteration &iter, r_sparseMat& A, std::vector<T> & x, const std::vector<T> & rhs, const std::vector<T> & xd, const std::vector<int>& ld);

/** conjugate gradient with diagonal preconditioner, returns residu. Template generic_cg is called with dummy xd and ld */
template <typename T>
T cg(iteration &iter, r_sparseMat& A, std::vector<T> & x, std::vector<T> & b)
    {
    std::vector<T> xd;
    std::vector<int> ld;
    return generic_cg<false,T>(iter,A,x,b,xd,ld);
    }

/** conjugate gradient with diagonal preconditioner with Dirichlet conditions, returns residu */
template <typename T>
T cg_dir(iteration &iter, r_sparseMat& A, std::vector<T> & x, const std::vector<T> & rhs,
         const std::vector<T> & xd, const std::vector<int>& ld)
    { return generic_cg<true,T>(iter,A,x,rhs,xd,ld); }

/** generic template for stabilized bigradient conjugate  */
template <bool MASK,typename T>
T generic_bicg(iteration &iter, r_sparseMat& A, std::vector<T> & x, const std::vector<T> & rhs,
               const std::vector<T>& xd, const std::vector<int>& ld);

/** biconjugate gradient and diagonal preconditionner, returns residu. Template generic_bicg is called with dummy xd and ld  */
template <typename T>
T bicg(iteration &iter, r_sparseMat& A, std::vector<T> & x, const std::vector<T> & b)
    {
    std::vector<T> xd;
    std::vector<int> ld;
    return generic_bicg<false,T>(iter,A,x,b,xd,ld);
    }

/** biconjugate gradient with diagonal preconditioner with Dirichlet conditions, returns residu */
template <typename T>
T bicg_dir(iteration &iter, r_sparseMat& A, std::vector<T> & x, const std::vector<T> & rhs,
           const std::vector<T> & xd, const std::vector<int>& ld)
    { return generic_bicg<true,T>(iter,A,x,rhs,xd,ld); }

/** operator<< for r_sparseVect */
inline std::ostream & operator<<(std::ostream & flux, r_sparseVect const& v)
    {v.print(flux); return flux;}

/** operator<< for r_sparseMat */
inline std::ostream & operator<<(std::ostream & flux, r_sparseMat const& m)
    {m.print(flux); return flux;}

} // end namespace algebra

#endif

