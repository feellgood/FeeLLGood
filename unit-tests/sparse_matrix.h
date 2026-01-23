/** \file sparse_matrix.h
 *
 * For testing purposes, we would like to build a sparse matrix in one step, from a list of
 * coefficients, instead of having to first build the matrix shape and then add the coefficients.
 * The function buildSparseMat() fulfills this purpose.
 */

#include <vector>

#include "../algebra/sparseMat.h"

using algebra::SparseMatrix;
using algebra::MatrixShape;

/** Matrix coefficient consumed by buildSparseMat() below. */
struct MatrixCoefficient
    {
    int row, col;
    double value;
    };

/** Helper function to build a one-off sparse matrix from a list of coefficients. */
inline SparseMatrix buildSparseMat(int N, const std::vector<MatrixCoefficient>& coefficients)
    {
    MatrixShape shape(N);
    for (const MatrixCoefficient& coefficient : coefficients)
        {
        shape[coefficient.row].insert(coefficient.col);
        }
    SparseMatrix matrix(shape);
    for (const MatrixCoefficient& coefficient : coefficients)
        {
        matrix.add(coefficient.row, coefficient.col, coefficient.value);
        }
    return matrix;
    }
