#ifndef solver_h
#define solver_h

/** \file solver.h
 \brief two templates to fill matrix and vectors in various dimensionnality situations. DIM_PROBLEM = 1 is used for electrostatics (V) DIM_PROBLEM = 3 is used for spin accumulation (Q has three components)
Warning : DIM_PROBLEM = 2 cannot be used for micromagnetic problem. The latter is solved in the tangent plane of the magnetization plane, leading to a different indices computation and matrix filling than here.
TODO: these templates could be specialized for DIM_PROBLEM = 2 (see warning above)
*/

#include <eigen3/Eigen/Dense>
#include "mesh.h"

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
                        const int max_iter /**< [in] maximum number of iterations */):
                        DIM_PB(DIM_PROBLEM), msh(&_msh), paramTet(_pTetra), paramFac(_pFac), iter(name,_tol,v,max_iter) {}

    protected:
        /** dimensionnality of th problem */
        const int DIM_PB;

        /** mesh pointer to access nodes, fac, tet, and others geometrical values and methods */
        Mesh::mesh *msh;

        /** this vector contains the material parameters for all volume regions for all the tetrahedrons */
        std::vector<Tetra::prm> paramTet;

        /** this vector contains the material parameters for all surface regions for all the triangular facettes */
        std::vector<Facette::prm> paramFac;

        /** monitor the solver called in method solve() */
        algebra::iteration<double> iter;

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
