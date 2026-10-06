#ifndef BLOCK_PRECOND_H
#define BLOCK_PRECOND_H

#include <stdexcept>
#include <vector>

#include <eigen3/Eigen/Dense>

#include "sparseMat.h"

namespace algebra
{
/** \class BlockDiagPrecond
block diagonal preconditioner: inverse of the BS x BS diagonal block of each node, for a system
whose unknowns are numbered node by node (unknown BS*i + d is the component d of the node i).
For the spin diffusion problem (BS = 3), this block holds the relaxation and the precession
coupling s x u, local to each node, which the diagonal preconditioner misses (same preconditioner
as the BLOCK3 of st-feeLLGood).
The rows and columns of the Dirichlet components are replaced by the identity before inversion, so
that the preconditioner does not mix the imposed components with the free ones.
*/
template <int BS>
class BlockDiagPrecond
    {
public:
    /** computes the inverse of the diagonal block of each node of A; ld: indices of the Dirichlet
     * components */
    BlockDiagPrecond(const SparseMatrix &A, const int nbBlocks, const std::vector<int> &ld)
        : inv(nbBlocks)
        {
        std::vector<char> dir(BS * nbBlocks, 0);
        for (const int i : ld)
            { dir[i] = 1; }
        for (int n = 0; n < nbBlocks; n++)
            {
            Eigen::Matrix<double, BS, BS> a;
            for (int r = 0; r < BS; r++)
                for (int c = 0; c < BS; c++)
                    {
                    const int I = BS * n + r, J = BS * n + c;
                    a(r, c) = (dir[I] || dir[J]) ? (r == c ? 1.0 : 0.0) : A(I, J);
                    }
            if (a.determinant() == 0.0)
                { throw std::runtime_error("BlockDiagPrecond: singular diagonal block"); }
            inv[n] = a.inverse();
            }
        }

    /** y = M x, M being the block diagonal preconditioner */
    void apply(const std::vector<double> &x, std::vector<double> &y) const
        {
        for (size_t n = 0; n < inv.size(); n++)
            {
            Eigen::Map<const Eigen::Matrix<double, BS, 1>> xb(&x[BS * n]);
            Eigen::Map<Eigen::Matrix<double, BS, 1>> yb(&y[BS * n]);
            yb = inv[n] * xb;
            }
        }

private:
    /** inverses of the diagonal blocks */
    std::vector<Eigen::Matrix<double, BS, BS>> inv;
    };

}  // namespace algebra
#endif
