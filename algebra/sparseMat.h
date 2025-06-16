#ifndef SPARSEMAT_H
#define SPARSEMAT_H

/** \file sparseMat.h
 \brief Sparse matrices

The main class provided here is `r_sparseMat`, a "read-mode" square sparse matrix. It is meant to
efficiently compute matrix-vector products, where the vectors are full (not sparse), represented as
std::vector<double>. Such a matrix can be constructed in two ways:

- If the set of valid index pairs (locations that can have non-zero values) is **not** known in
  advance, one can construct a `w_sparseMat` (write-mode sparse matrix), populate it with its
  `insert()` method, then use it to construct the read-mode sparse matrix.

- If the set of valid index pairs is known, one can store it in a `MatrixShape`, then use this shape
  to construct the read-mode sparse matrix with all values initialized to zero.

Once constructed, the set of valid index pairs is immutable. However, the associated values can be
modified with `clear()` and `add()`. Rewriting the values this way is way more efficient than
reconstructing the matrix from scratch.
 */

#include <set>
#include <map>
#include <vector>
#include <iostream>
#include <cassert>
#include <algorithm>
#include <execution>

#include "algebraCore.h"

namespace algebra
{

/**
The shape of a sparse matrix is the set of valid indices. Specifically, the shape element at index i
is the set of j indices such that (i, j) is a valid index pair for the matrix.
 */
using MatrixShape = std::vector<std::set<int>>;

/**
\class w_sparseMat
\brief Write-mode sparse matrix.

This is a container for matrix elements that can be modified by inserting elements at arbitrary
locations.
*/
class w_sparseMat
{
    friend class r_sparseMat;

    /**
    Write-mode sparse vector. This is an (index -> value) map.
    */
    using w_sparseVect = std::map<int, double>;

public:
    /** Constructor. N is the matrix size: this is a square N×N matrix. */
    w_sparseMat(int N) : m(N) {}

    /**
    Insert a matrix element with value val, at line i, column j. If an element already exists at
    this position, add `val` to its value.

    The caller must ensure both indices lie within [0, N).
    */
    void insert(const int i, const int j, const double val)
        {
        assert(i >= 0 && i < m.size());
        assert(j >= 0 && j < m.size());
        m[i][j] += val;
        }

private:
    std::vector<w_sparseVect> m;
};


/** \class r_sparseMat
\brief Read-mode square sparse matrix.

This is the main sparse matrix class, meant for efficiently computing matrix-vector products.
*/
class r_sparseMat
{
    /**
    Read-mode sparse vector.

    This holds a list of indices, and a list of matching values.
     - both lists have the same size
     - the indices are unique and sorted
     - the values match the indices in the sense that `values[k]` is the vector element at index
       `indices[k]`.
    */
    struct r_sparseVect
    {
        std::vector<int> indices;  /**< array of vector indices. */
        std::vector<double> values;  /**< array of vector values matching the indices */
    };

public:
    /** Construct from a read-mode sparse matrix. */
    r_sparseMat(const w_sparseMat &A)
        {
        m.reserve(A.m.size());
        for (const w_sparseMat::w_sparseVect& line_in: A.m)
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

    /** Construct from a matrix shape. Initialize all stored values to zero. */
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

    /** Zero-out all elements, while preserving the shape. */
    void clear()
        {
        for (r_sparseVect& line: m)
            for (double& value: line.values)
                value = 0;
        }

    /** Add the given value at position (i, j), which must belong to the shape. */
    void add(int i, int j, double val)
        {
        assert(i >= 0 && i < m.size());
        assert(j >= 0 && j < m.size());
        r_sparseVect& line = m[i];
        auto it = std::lower_bound(line.indices.begin(), line.indices.end(), j);
        assert(it != line.indices.end() && *it == j);
        int k = it - line.indices.begin();
        line.values[k] += val;
        }

    /** Print to the provided stream. */
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

    /** Return the value at position (i, j). */
    double operator() (const int i, const int j) const
        {
        assert(i >= 0 && i < m.size());
        assert(j >= 0 && j < m.size());
        const r_sparseVect& line = m[i];
        auto it = std::lower_bound(line.indices.begin(), line.indices.end(), j);
        if (it == line.indices.end() || *it != j)
            return 0;
        int k = it - line.indices.begin();
        return line.values[k];
        }

    /** Y = this*X */
    template <typename T>
    void mult(std::vector<T> const& X, std::vector<T> &Y)
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

    /** Build a diagonal preconditioner D from this matrix. */
    template <typename T>
    void build_diag_precond(std::vector<T> &D) const
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
