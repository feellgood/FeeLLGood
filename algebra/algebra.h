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
/** Y *= alpha */
inline void scaled( const double alpha, std::vector<double> & Y) 
	{ std::for_each(Y.begin(),Y.end(),[alpha](double &_x){ _x *= alpha; }); }

/** returns scalar product X.Y */
inline double dot(const std::vector<double> & X,const std::vector<double> & Y)
	{ return std::inner_product(X.begin(),X.end(),Y.begin(),0.0); }

/** direct product : Z = X⊗Y */
inline void p_direct(const std::vector<double> & X,const std::vector<double> & Y,std::vector<double> & Z)
	{ for(unsigned int i=0;i<Z.size();i++) Z[i]=X[i]*Y[i]; }

/** Y += X       */
inline void add(const std::vector<double> & X, std::vector<double> & Y)
	{ std::transform(Y.begin(),Y.end(),X.begin(),Y.begin(),std::plus<double>()  ); }

/** Y -= X       */
inline void sub(const std::vector<double> & X, std::vector<double> & Y)
	{ std::transform(Y.begin(),Y.end(),X.begin(),Y.begin(),std::minus<double>()  ); }

/** Y += alpha*X       */
inline void scaled_add(const std::vector<double> & X,const double alpha, std::vector<double> & Y)
	{ std::transform(Y.begin(),Y.end(),X.begin(),Y.begin(),[alpha] (const double _x,double _y) { return _x+(alpha*_y); }   ); }

/** euclidian norm of vector X */
inline double norm(const std::vector<double> & X)
	{ return sqrt(fabs( dot(X,X) )); }

/** Y = A*X with r_sparseMat A */
inline void mult(r_sparseMat & A,std::vector<double> const& X,std::vector<double> &Y)
    {
    const int _size = X.size();
    Y.resize(_size);
    if (A.getDim() == _size)
	    { for(int i=0;i<_size;i++) Y[i]= A(i).dot(X); }
    }

/** apply a mask to vector X : all coefficients in vector mask are zeroed in X */
inline void applyMask(const std::vector<int>& mask, std::vector<double> & X)
    { std::for_each(mask.begin(),mask.end(),[&X](const int _i){ X[_i] = 0.0; }); }

/** generic template for gradient conjugate */
template<bool MASK>
double generic_cg(iteration &iter, r_sparseMat& A, std::vector<double> & x, const std::vector<double> & rhs, const std::vector<double> & xd, const std::vector<int>& ld);

/** conjugate gradient with diagonal preconditioner, returns residu. Template generic_cg is called with dummy xd and ld */
inline double cg(r_sparseMat& A, std::vector<double> & x, const std::vector<double> & b, iteration &iter)
    {
    std::vector<double> xd;
    std::vector<int> ld;
    return generic_cg<false>(iter,A,x,b,xd,ld);
    }

/** conjugate gradient with diagonal preconditioner with Dirichlet conditions, returns residu */
inline double cg_dir(r_sparseMat& A, std::vector<double> & x, const std::vector<double> & rhs, 
              const std::vector<double> & xd, const std::vector<int>& ld, iteration &iter)
    { return generic_cg<true>(iter,A,x,rhs,xd,ld); }

/** generic template for stabilized bigradient conjugate  */
template<bool MASK>
double generic_bicg( iteration &iter, r_sparseMat& A, std::vector<double> & x, const std::vector<double> & rhs, const std::vector<double>& xd, const std::vector<int>& ld);

/** biconjugate gradient and diagonal preconditionner, returns residu. Template generic_bicg is called with dummy xd and ld  */
inline double bicg(r_sparseMat& A, std::vector<double> & x, const std::vector<double> & b, iteration &iter)
    {
    std::vector<double> xd;
    std::vector<int> ld;
    return generic_bicg<false>(iter,A,x,b,xd,ld);
    }

/** biconjugate gradient with diagonal preconditioner with Dirichlet conditions, returns residu */
inline double bicg_dir(r_sparseMat& A, std::vector<double> & x, const std::vector<double> & rhs, 
              const std::vector<double> & xd, const std::vector<int>& ld, iteration &iter)
    { return generic_bicg<true>(iter,A,x,rhs,xd,ld); }

/** operator<< for r_sparseVect */
inline std::ostream & operator<<(std::ostream & flux, r_sparseVect const& v)
    {v.print(flux); return flux;}

/** operator<< for r_sparseMat */
inline std::ostream & operator<<(std::ostream & flux, r_sparseMat const& m)
    {m.print(flux); return flux;}

} // end namespace algebra

#endif

