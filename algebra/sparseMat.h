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
#include <mutex>

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
    using SparseVector = std::map<int, double>;

public:
    /** Constructor. N is the matrix size: this is a square N×N matrix. */
    w_sparseMat(int N) : rows(N) {}

    /**
    Insert a matrix element with value val, at row i, column j. If an element already exists at this
    position, add `val` to its value.

    The caller must ensure both indices lie within [0, N).
    */
    void insert(const int i, const int j, const double val)
        {
        assert(i >= 0 && i < rows.size());
        assert(j >= 0 && j < rows.size());
        rows[i][j] += val;
        }

private:
    std::vector<SparseVector> rows;
};


/** \class r_sparseMat
\brief Read-mode square sparse matrix.

This is the main sparse matrix class, meant for efficiently computing matrix-vector products.
*/
class r_sparseMat
{
    /**
    Read-mode sparse vector.

    This holds a mutex that protects the data during calls to add(), a list of indices, and a list
    of matching values. The following invariants should be enforced:
     - both lists have the same size
     - the indices are unique and sorted
     - the values match the indices in the sense that `values[k]` is the vector element at index
       `indices[k]`.

    Since mutexes are not movable, neither is this class.
    */
    struct SparseVector
    {
        std::mutex mutex;   /**< mutex guarding the values during calls to add() */
        std::vector<int> indices;  /**< array of vector indices. */
        std::vector<double> values;  /**< array of vector values matching the indices */
    };

public:
    /** Construct from a read-mode sparse matrix. */
    r_sparseMat(const w_sparseMat &source) : rows(source.rows.size())
        {
        for (size_t i = 0; i < source.rows.size(); ++i)
            {
            const w_sparseMat::SparseVector& source_row = source.rows[i];
            SparseVector& row = rows[i];
            row.indices.reserve(source_row.size());
            row.values.reserve(source_row.size());
            for (auto it = source_row.begin(); it != source_row.end(); ++it)
                {
                row.indices.push_back(it->first);
                row.values.push_back(it->second);
                }
            }
        }

    /** Construct from a matrix shape. Initialize all stored values to zero. */
    r_sparseMat(const MatrixShape &shape) : rows(shape.size())
        {
        for (size_t i = 0; i < shape.size(); ++i)
            {
            const std::set<int>& row_shape = shape[i];
            SparseVector& row = rows[i];
            row.indices.reserve(row_shape.size());
            for (auto it = row_shape.begin(); it != row_shape.end(); ++it)
                row.indices.push_back(*it);
            row.values.resize(row_shape.size());
            }
        }

    /** Zero-out all elements, while preserving the shape. */
    void clear()
        {
        for (SparseVector& row: rows)
            for (double& value: row.values)
                value = 0;
        }

    /**
    Add the provided value at position (i, j), which must belong to the shape. This method is
    thread-safe.
    */
    void add(int i, int j, double val)
        {
        assert(i >= 0 && i < rows.size());
        assert(j >= 0 && j < rows.size());
        SparseVector& row = rows[i];
        auto it = std::lower_bound(row.indices.begin(), row.indices.end(), j);
        assert(it != row.indices.end() && *it == j);
        int k = it - row.indices.begin();
        std::mutex& mutex = row.mutex;
        double& value = row.values[k];
        mutex.lock();
        value += val;
        mutex.unlock();
        }

    /** Print to the provided stream. */
    void print(std::ostream & flux = std::cout) const
        {
        flux << "[\n";
        for (const SparseVector& row: rows)
            {
            flux << "  {";
            for (size_t k = 0; k < row.indices.size(); ++k)
                {
                flux << row.indices[k] << ": " << row.values[k];
                flux << (k < row.indices.size()-1 ? ", " : "\n");
                }
            flux << "}\n";
            }
        flux << ']';
        }

    /** Return the value at position (i, j). */
    double operator() (const int i, const int j) const
        {
        assert(i >= 0 && i < rows.size());
        assert(j >= 0 && j < rows.size());
        const SparseVector& row = rows[i];
        auto it = std::lower_bound(row.indices.begin(), row.indices.end(), j);
        if (it == row.indices.end() || *it != j)
            return 0;
        int k = it - row.indices.begin();
        return row.values[k];
        }

    /** Y = this*X */
    template <typename T>
    void mult(std::vector<T> const& X, std::vector<T> &Y)
        {
        assert(X.size() == rows.size());
        assert(Y.size() == rows.size());
        std::transform(EXEC_POL, rows.begin(), rows.end(), Y.begin(),
            [&X](const SparseVector &row)
            {
            double val = 0;
            for (size_t k = 0; k < row.indices.size(); ++k)
                val += row.values[k] * X[row.indices[k]];
            return val;
            });
        }

    /** Build a diagonal preconditioner D from this matrix. */
    template <typename T>
    void build_diag_precond(std::vector<T> &D) const
        {
        assert(D.size() == rows.size());
        for (size_t i = 0; i < rows.size(); ++i)
            {
            const double c = (*this)(i,i);
            assert(c != 0);
            D[i] = 1.0 / c;
            }
        }

private:
    /**
    Container for the matrix rows. Since SparseVector is not movable, this container cannot be
    resized.
    */
    std::vector<SparseVector> rows;
}; // end class r_sparseMat

} // end namespace algebra

#endif
