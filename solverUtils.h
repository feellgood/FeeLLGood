#ifndef solverUtils_h
#define solverUtils_h

/** \file solverUtils.h
 \brief two templates to fill matrix and vectors in various dimensionnality situations. DIM_PROBLEM = 1 is used for electrostatics (V) DIM_PROBLEM = 3 is used for spin accumulation (Q has three components)
Warning : DIM_PROBLEM = 2 cannot be used for micromagnetic problem. The latter is solved in the tangent plane of the magnetization plane, leading to a different indices computation and matrix filling than here.
TODO: these templates could be specialized for DIM_PROBLEM = 2 (see warning above)
*/

#include <eigen3/Eigen/Dense>

/** function template.
first parameter N is the number of indices of the element to build matrix from: ind.size() = N
second parameter DIM_PROBLEM is the dimension of the quantity to solve (should be 1,2 or 3)  */
template <int N, int DIM_PROBLEM>
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
                    //K.insert(di*NOD+i_, dj*NOD+j_, Ke(di*N+ie,dj*N+je));
            }
        }
    }

/** function template.
first parameter N is the number of indices of the element to build vector from
second parameter DIM_PROBLEM is the dimension of the quantity to solve (should be 1,2 or 3)  */
template <int N, int DIM_PROBLEM>
void buildVect(std::vector<int> &ind, std::vector<double> &Le, std::vector<double> &L)
    {
    for (int ie=0; ie<N; ie++)
        {
        int i_ = ind[ie];
        for (int di=0; di<DIM_PROBLEM; di++)
            { L[DIM_PROBLEM*i_ + di] += Le[di*N+ie]; }//{ L[di*NOD+i_] += Le[di*N+ie]; }
        }
    }

/** template class to describe boundary conditions
 * T should be either a 3D vector or a double
 * indices stored in facIdx are the indices defining a surface in mesh facette list where the
 * boundary condition does apply
 * */
template <class T>
class boundaryCondition
    {
    public:
    std::vector<int> facIdx; /**< list of the facette indices defining a surface stored in mesh */
    T value; /**< value on the surface defined by the list facIdx */
    };

#endif
