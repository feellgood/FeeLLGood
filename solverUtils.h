/** \file solverUtils.h
 \brief two templates to fill matrix and vectors in various dimensionnality situations. DIM_PROBLEM = 1 is used for electrostatics (V) DIM_PROBLEM = 3 is used for spin accumulation (Q has three components)
Warning : DIM_PROBLEM = 2 cannot be used for micromagnetic problem. The latter is solved in the tangent plane of the magnetization plane, leading to a different indices computation and matrix filling than here.
TODO: these templates could be specialized for DIM_PROBLEM = 2 (see warning above)
*/

/** function template.
first parameter N is the number of indices of the element to build matrix from: ind.size() = N
second parameter DIM_PROBLEM is the dimension of the quantity to solve (should be 1,2 or 3)  */
template <int N, int DIM_PROBLEM>
void buildMat(const int NOD, std::vector<int> &ind,
              Eigen::Matrix<double,DIM_PROBLEM*N,DIM_PROBLEM*N> &Ke, algebra::w_sparseMat &K)
    {
    for (int ie=0; ie<N; ie++)
        {
        int i= ind[ie];
        for (int je=0; je<N; je++)
            {
            int j= ind[je];
            for (int di=0; di<DIM_PROBLEM; di++)
                for (int dj=0; dj<DIM_PROBLEM; dj++)
                    K.insert(di*NOD+i, dj*NOD+j, Ke(di*N+ie,dj*N+je));
            }
        }
    }

/** function template.
first parameter N is the number of indices of the element to build vector from
second parameter DIM_PROBLEM is the dimension of the quantity to solve (should be 1,2 or 3)  */
template <int N, int DIM_PROBLEM>
void buildVect(const int NOD, std::vector<int> &ind,
               std::vector<double> &Le, std::vector<double> &L)
    {
    for (int ie=0; ie<N; ie++)
        {
        int i= ind[ie];
        for (int di=0; di<DIM_PROBLEM; di++) { L[di*NOD+i] += Le[di*N+ie]; }
        }
    }
