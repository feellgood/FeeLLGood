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
class w_sparseMat : private std::vector<w_sparseVect>
{
    friend class r_sparseMat;

public:
    /** constructor */
    w_sparseMat(int N) : std::vector<w_sparseVect>(N) {}

    /** inserter for a coefficient val at line i col j */
    void insert(const int i, const int j, const double val)
        {
        assert(i >= 0 && i < size());
        assert(j >= 0 && j < size());
        (*this)[i][j] += val;
        }

}; // end class w_sparseMat


/** \class r_sparseMat
read sparse matrix
The constructor is building the data from a write sparse matrix to access efficiently the coefficients values
*/
class r_sparseMat
{
public:
    /** construct a sparse matrix from a "read mode" matrix */
    r_sparseMat(const w_sparseMat &A)
        {
        m.reserve(A.size());  // A.size() is the number of lines
        for (const w_sparseVect& line_in: A)
            {
            r_sparseVect line;
            line.indices.reserve(line_in.size());
            line.values.reserve(line_in.size());
            for (auto it = line_in.begin(); it != line_in.end(); ++it)
                {
                line.indices.push_back(it->first);
                line.values.push_back(it->second);
                }
            m.push_back(line);
            }
        }

    /** construct a sparse matrix from its shape */
    r_sparseMat(const MatrixShape &shape)
        {
        m.reserve(shape.size());
        for (const std::set<int>& line_shape: shape)
            {
            r_sparseVect line;
            line.indices.reserve(line_shape.size());
            for (auto it = line_shape.begin(); it != line_shape.end(); ++it)
                line.indices.push_back(*it);
            line.values.resize(line_shape.size());
            m.push_back(line);
            }
        }

    /** zero all elements, while preserving the shape */
    void clear()
        {
        for (r_sparseVect& line: m)
            for (double& value: line.values)
                value = 0;
        }

    /** add the value at position (i, j), which must belog to the shape */
    void add(int i, int j, double val)
        {
        assert(i >= 0 && i < m.size());
        assert(j >= 0 && j < m.size());
        r_sparseVect& line = m[i];
        size_t k = 0;
        for (; k < line.indices.size(); ++k)
            if (line.indices[k] == j)
                break;
        assert(k < line.indices.size());
        line.values[k] += val;
        }

    /** printing function */
    void print(std::ostream & flux = std::cout) const
        {
        flux << "[\n";
        for (const r_sparseVect& line: m)
            {
            flux << "  {";
            for (size_t k = 0; k < line.indices.size(); ++k)
                {
                flux << line.indices[k] << ": " << line.values[k];
                flux << (k < line.indices.size()-1 ? ", " : "\n");
                }
            flux << "}\n";
            }
        flux << ']';
        }

    /** getter for a coefficient value */
    double operator() (const int i, const int j) const
        {
        assert(i >= 0 && i < m.size());
        assert(j >= 0 && j < m.size());
        const r_sparseVect& line = m[i];
        for (size_t k = 0; k < line.indices.size(); ++k)
            if (line.indices[k] == j)
                return line.values[k];
        return 0;
        }

    /** Y = this*X */
    template <typename T>
    void mult(Vector<T> const& X, Vector<T> &Y)
        {
        assert(X.size() == m.size());
        assert(Y.size() == m.size());
        std::transform(std::execution::par, m.begin(), m.end(), Y.begin(),
            [&X](const r_sparseVect &line)
            {
            double val = 0;
            for (size_t k = 0; k < line.indices.size(); ++k)
                val += line.values[k] * X[line.indices[k]];
            return val;
            });
        }

    /** build diagonal preconditioner D from input matrix(this) */
    template <typename T>
    void build_diag_precond(Vector<T> &D) const
        {
        assert(D.size() == m.size());
        for (size_t i = 0; i < m.size(); ++i)
            {
            const double c = (*this)(i,i);
            assert(c != 0);
            D[i] = 1.0 / c;
            }
        }

private:
    /** coefficient container */
    std::vector<r_sparseVect> m;
}; // end class r_sparseMat



} // end namespace algebra

#endif
