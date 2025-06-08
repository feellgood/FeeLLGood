#ifndef SPARSEMAT_H
#define SPARSEMAT_H

/** \file sparseMat.h
 \brief read and write sparse matrix
r_sparseMat : read sparse matrix : it is buit calling the constructor with a w_sparseMat as argument
w_sparseMat : write sparse matrix : a std::vector of w_sparse_vector of coefficients (i,value)
 */

#include <iostream>
#include <set>
#include <vector>
#include <cassert>
#include <execution>

#include "algebraCore.h"

namespace algebra
{

/**
The shape of a sparse matrix is the set of valid indices.
 */
using MatrixShape = std::vector<std::set<int>>;

/**
\class w_sparseMat
write sparse Matrix, it is a container for coefficients of a 'line' sparse matrix.
If some m_coeff have the same indices, they will be summed to build the final matrix coefficient
*/
class w_sparseMat
{
    friend class r_sparseMat;

public:
/** constructor */
    w_sparseMat(const int _N) : N(_N)
        {
        m.reserve(N);
        for (int i = 0; i < N; ++i)
            m.emplace_back(N);
        }

/** inserter for a coefficient val at line i col j */
    void insert(const int i, const int j, const double val)
        {
        assert(i<N && j<N );
        m[i].insert(j, val);
        }

/** getter for the number of lines */
    int getDim() const {return N;}

private:
/** dimension of sparse matrix, N is the number of lines */
    const int N;

/** container for the write sparse matrix coefficients */
    std::vector<w_sparseVect> m;
}; // end class w_sparseMat


/** \class r_sparseMat
read sparse matrix
The constructor is building the data from a write sparse matrix to access efficiently the coefficients values
*/
class r_sparseMat
{
public:
    /** construct a sparse matrix from a "read mode" matrix */
    r_sparseMat(const w_sparseMat &A):N(A.getDim())
        {
        m.reserve(N);  // N is the number of lines
        if (!A.m.empty())
            { for(int i=0; i<N; ++i){ m.emplace_back(A.m[i]); } }
        }

    /** construct a sparse matrix from its shape */
    r_sparseMat(const MatrixShape &shape) : N(shape.size())
        {
        for (int i = 0; i < N; ++i)
            { m.emplace_back(N, shape[i]); }
        }

    /** zero all elements, while preserving the shape */
    void clear()
        {
        for (r_sparseVect& line: m)
            line.clear();
        }

    /** add the value at position (i, j), which must belog to the shape */
    void add(int i, int j, double val)
        {
        m[i].add(j, val);
        }

/** printing function */
    void print(std::ostream & flux = std::cout) const
    { std::for_each(m.begin(),m.end(),[&flux](r_sparseVect const& _v) {_v.print(flux);} ); }

    /** getter for the number of lines */
    int getDim() const {return N;}

/** getter for an innner sparse vector */
    const r_sparseVect & operator() (const int i) const {return m[i];}

/** getter for a coefficient value */
    double operator() (const int i, const int j) const { return m[i].getVal(j); }

/** Y = this*X */
template <typename T>
void mult(Vector<T> const& X, Vector<T> &Y)
    {
    if (int(X.size()) != N)
        { std::cerr << "Error: wrong dimensions in matrix-vector product.\n"; exit(1); }
    std::transform(std::execution::par,m.begin(),m.end(),Y.begin(),
                   [&X](const r_sparseVect &_v){ return _v.dot(X); });
    }

/** build diagonal preconditioner D from input matrix(this) */
template <typename T>
void build_diag_precond(Vector<T> &D) const
    {
    for(unsigned int i=0;i<D.size();i++)
        {
        const double c = (*this)(i,i);
        if(c != 0)
            { D[i] = 1.0/c; }
        else
            { std::cout <<"Error: zero on sparse matrix diagonal.\n"; exit(1); }
        }
    }

private:
/** dimension of the sparse matrix (nb of lines) */
    const int N;

/** coefficient container */
    std::vector<r_sparseVect> m;
}; // end class r_sparseMat



} // end namespace algebra

#endif
