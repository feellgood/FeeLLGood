#ifndef ALG_H
#define ALG_H

/** \file alg.h 
 * \brief set of class to handle sparse matrix operations for gradient conjugate algorithm
 * a sparse vector class
 * a read and a write sparse matrix class
 * various vector operations, scalar and direct products; ...
 * most of the member functions of the classes are using lambdas and C++11 algorithm and numeric
 * */

#include <iostream>
#include <cmath> // sqrt,fabs
#include <algorithm>
#include <numeric> // inner_product

#include "coeff.h"
#include "sparseVect.h"
#include "sparseMat.h"
#include "denseMat.h"
#include "iter.h"

/** \namespace alg
 * grab altogether sparse matrix, vector and dedicated functions
*/

namespace alg
{
/** 
Y = alpha*X 
*/
inline void scaled(const std::vector<double> & X, const double alpha, std::vector<double> & Y) 
	{ std::transform(X.begin(),X.end(),Y.begin(),[alpha](const double _x){ return alpha*_x; }); }

/** 
Y *= alpha 
*/
inline void scaled( const double alpha, std::vector<double> & Y) 
	{ std::for_each(Y.begin(),Y.end(),[alpha](double &_x){ _x *= alpha; }); }

/**
returns scalar product X.Y
*/
inline double dot(const std::vector<double> & X,const std::vector<double> & Y)
	{ return std::inner_product(X.begin(),X.end(),Y.begin(),0.0); }

/**
direct product : Z = X⊗Y
*/
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
	{ return sqrt(fabs( alg::dot(X,X) )); }

/** Y = A*X with sparseMat A */
inline void mult(alg::r_sparseMat & A,std::vector<double> const& X,std::vector<double> &Y)
    {
    const int _size = X.size();
    Y.resize(_size);
    if (A.getDim() == _size)
	    { for(int i=0;i<_size;i++) Y[i]= A(i).dot(X); }
    }

/** Y = A*X with denseMat A */
inline void mult(alg::denseMat const& A,std::vector<double> const& X,std::vector<double> &Y)
{
const size_t _size = X.size();
assert(A.ncols() == _size);

Y.resize(_size);
for (size_t i=0; i<A.nrows(); i++) {
       double val(0);
       for (size_t j=0; j<_size; j++) val += A(i,j)*X[j];
       Y[i]=val;
       }
}

/** Y = trans(X)*A */
inline void transposed_mult(std::vector<double> const& X,alg::denseMat const& A,std::vector<double> &Y)
{
const size_t ncolsA = A.ncols();
Y.resize(ncolsA);
const size_t _size = X.size();

for (size_t j=0; j<ncolsA; j++) { 
       double val(0);
       for (size_t i=0; i<_size; i++) val += X[i]*A(i,j);
       Y[j]=val;
       }
}

/** Y = trans(A)*X with denseMat A */
inline void transposed_mult(alg::denseMat const& A,std::vector<double> const& X,std::vector<double> &Y)
{
const size_t nrowsA = A.nrows();
assert(nrowsA == X.size());

Y.resize(A.ncols());
for (size_t i=0; i<A.ncols(); i++) {
       double val(0);
       for (size_t j=0; j<nrowsA; j++) val += A(j, i)*X[j];
       Y[i]=val;
       }
}


/** C = trans(A)*B */
inline void transposed_mult(alg::denseMat const& A,alg::denseMat const& B,alg::denseMat & C)
{
const size_t nrowsA = A.nrows();
const size_t ncolsA = A.ncols();
const size_t ncolsB = B.ncols();
assert(nrowsA==B.nrows());

C.nrows()=ncolsA;
C.ncols()=ncolsB;
C.resize(ncolsA*ncolsB);
for (size_t i=0; i<ncolsA; i++) { 
    for (size_t j=0; j<ncolsB; j++) 
	{
        double val(0);
        for (size_t k=0; k<nrowsA; k++) { val += A(k, i)*B(k, j); }
	C(i, j) = val;      
	}
}
}

inline void applyMask(const std::vector<int>& mask, std::vector<double> & x)
    { std::for_each(mask.begin(),mask.end(),[&x](const int _i){ x[_i] = 0.0; }); }

template<bool MASK>
double generic_cg(alg::iteration &iter, alg::r_sparseMat& A, std::vector<double> & x, const std::vector<double> & rhs, const std::vector<double> & xd, const std::vector<int>& ld);

/** conjugate gradient with diagonal preconditioner, returns residu */
inline double cg(alg::r_sparseMat& A, std::vector<double> & x, const std::vector<double> & b, alg::iteration &iter)
    {
    std::vector<double> xd;
    std::vector<int> ld;
    return alg::generic_cg<false>(iter,A,x,b,xd,ld);
    }

/** conjugate gradient with diagonal preconditioner with Dirichlet conditions, returns residu */
inline double cg_dir(alg::r_sparseMat& A, std::vector<double> & x, const std::vector<double> & rhs, 
              const std::vector<double> & xd, const std::vector<int>& ld, alg::iteration &iter)
    {
    return alg::generic_cg<true>(iter,A,x,rhs,xd,ld);
    }

template<bool MASK>
double generic_bicg( alg::iteration &iter, alg::r_sparseMat& A, std::vector<double> & x, const std::vector<double> & rhs, const std::vector<double>& xd, const std::vector<int>& ld);

/** biconjugate gradient and diagonal preconditionner, returns residu */
inline double bicg(alg::r_sparseMat& A, std::vector<double> & x, const std::vector<double> & b, alg::iteration &iter)
    {
    std::vector<double> xd;
    std::vector<int> ld;
    return alg::generic_bicg<false>(iter,A,x,b,xd,ld);
    }

/** biconjugate gradient with diagonal preconditioner with Dirichlet conditions, returns residu */
inline double bicg_dir(alg::r_sparseMat& A, std::vector<double> & x, const std::vector<double> & rhs, 
              const std::vector<double> & xd, const std::vector<int>& ld, alg::iteration &iter)
    {
    return alg::generic_bicg<true>(iter,A,x,rhs,xd,ld);
    }
} // end namespace alg

#endif //ALG_H

