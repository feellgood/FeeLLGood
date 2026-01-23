#ifndef solver_h
#define solver_h

/** \file solver.h
 \brief two templates to fill matrix and vectors in various dimensionnality situations. DIM_PROBLEM = 1 is used for electrostatics (V) DIM_PROBLEM = 3 is used for spin accumulation (Q has three components)
Warning : DIM_PROBLEM = 2 cannot be used for micromagnetic problem. The latter is solved in the tangent plane of the magnetization plane, leading to a different indices computation and matrix filling than here.
TODO: these templates could be specialized for DIM_PROBLEM = 2 (see warning above)
*/

#include <eigen3/Eigen/Dense>
#include "mesh.h"
#include "algebra/algebra.h"

/** \class solver
 * \brief template class for the different solvers
 * template parameter DIM_PROBLEM: dimensionnality of the problem to solve
 * */

template <int DIM_PROBLEM>
class solver
    {
    public:
        /** constructor */
        explicit solver(Mesh::mesh & _msh /**< [in] mesh */,
                        std::vector<Tetra::prm> & _pTetra /**< [in] ref to vector of param tet (volume region parameters) */,
                        std::vector<Facette::prm> & _pFac /**< [in] ref to vector of param fac (surface region parameters) */,
                        const std::string name /**< [in] name of the solver method */,
                        const double _tol /**< [in] solver tolerance */,
                        const bool v /**< [in] verbose mode for iteration monitor */,
                        const int max_iter /**< [in] maximum number of iterations */,
                        std::function<bool(Mesh::Edge)> edge_filter = [](Mesh::Edge){ return true; }
                            /** [in] predicate for relevant mesh edges */):
                        msh(&_msh), NOD(_msh.getNbNodes()), paramTet(_pTetra),
                        paramFac(_pFac), verbose(v), iter(name,_tol,v,max_iter),
                        K(build_shape(edge_filter)), L_rhs(DIM_PROBLEM*NOD) {}

        /** check boundary conditions, exit if there is a mistake in the boundary conditions */
        virtual void checkBoundaryConditions(void) const = 0;

    protected:
        /** dimensionnality of the problem */
        static const int DIM_PB = DIM_PROBLEM;

        /** mesh pointer to access nodes, fac, tet, and others geometrical values and methods */
        Mesh::mesh *msh;

        /** number of nodes in the mesh */
        const int NOD;

        /** this vector contains the material parameters for all volume regions for all the tetrahedrons */
        const std::vector<Tetra::prm> &paramTet;

        /** this vector contains the material parameters for all surface regions for all the triangular facettes */
        const std::vector<Facette::prm> &paramFac;

        /** if verbose set to true, some printing are sent to terminal */
        const bool verbose;

        /** monitor the solver called in method solve() */
        algebra::iteration<double> iter;

        /** matrix of the system to solve */
        algebra::r_sparseMat K;

        /** RHS vector of the system to solve */
        std::vector<double> L_rhs;

        /** Build a matrix shape suitable for the current problem. */
        algebra::MatrixShape build_shape(std::function<bool(Mesh::Edge)> edge_filter)
            {
            algebra::MatrixShape shape(DIM_PROBLEM * NOD);

            // Add a DIM_PROBLEM × DIM_PROBLEM block connecting nodes i and j.
            auto add_block = [this, &shape](int i, int j)
                {
                for (int k = 0; k < DIM_PROBLEM; ++k)
                    {
                    for (int l = 0; l < DIM_PROBLEM; ++l)
                        { shape[DIM_PROBLEM*i+k].insert(DIM_PROBLEM*j+l); }
                    }
                };

            // Add a diagonal block for each node.
            for (int i = 0; i < NOD; ++i)
                { add_block(i, i); }

            // Add two off-diagonal blocks for each edge relevant to the current problem.
            for (auto edge: msh->edges)
                {
                if (edge_filter(edge))
                    {
                    add_block(edge.first, edge.second);
                    add_block(edge.second, edge.first);
                    }
                }

            return shape;
            }

        /** function template.
        parameter N is the number of indices of the element to build matrix from: ind.size() = N
         */
        template <int N>
        void buildMat(std::vector<int> &ind, Eigen::Matrix<double,DIM_PROBLEM*N,DIM_PROBLEM*N> &Ke, algebra::w_sparseMat &K)
            {
            for (int ie=0; ie<N; ie++)
                {
                int i_ = ind[ie];
                for (int je=0; je<N; je++)
                    {
                    int j_ = ind[je];
                    for (int di=0; di<DIM_PROBLEM; di++)
                        for (int dj=0; dj<DIM_PROBLEM; dj++)
                            K.insert(DIM_PROBLEM*i_ + di, DIM_PROBLEM*j_ + dj, Ke(di*N+ie,dj*N+je));
                    }
                }
            }

        /** function template.
        parameter N is the number of indices of the element to build vector from
        */
        template <int N>
        void buildVect(std::vector<int> &ind, std::vector<double> &Le, std::vector<double> &L)
            {
            for (int ie=0; ie<N; ie++)
                {
                int i_ = ind[ie];
                for (int di=0; di<DIM_PROBLEM; di++)
                    { L[DIM_PROBLEM*i_ + di] += Le[di*N+ie]; }
                }
            }
    }; // end template class solver

#endif
