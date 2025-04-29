#ifndef spinAccumulationSolver_h
#define spinAccumulationSolver_h

#include <vector>
#include "config.h"
#include "node.h"
#include "tetra.h"
#include "electrostatSolver.h"

/** \class spinAcc
 container for Spin Accumulation constants, diffusive accumulation spin model, Boundary conditions
 (potential fixed value on one surface, current density on another surface)
 */

// _tol could be 1e-6
class spinAcc
    {
    /** constructor */
    spinAcc(Mesh::mesh &_msh /**< [in] ref to the mesh */,
    electrostatSolver &_elec /**< [in] ref to the the electrostatic sub_problem */,
    std::vector<Tetra::prm> _pTetra /**< [in] ref to vector of param tetra (volume region parameters) */,
    std::vector<Facette::prm> _pFac /**< [in] ref to vector of param facette (surface region parameters) */,
    const double _tol /**< [in] tolerance for solvers */,
    const bool v /**< [in] verbose bool */,
    const int max_iter /**< [in] maximum number of iteration */):
        msh(_msh), elec(_elec), paramTetra(_pTetra),  paramFacette(_pFac), verbose(v), MAXITER(max_iter), NOD(_msh.getNbNodes())
        {
        Qs.resize(NOD);
        bool has_converged = solve(_tol);
        if (!has_converged)
            { std::cout << "spin accumulation solver has not converged." << std::endl; exit(1); }
        }

    private:
    /** mesh object to store nodes, fac, tet, and others geometrical values related to the mesh */
    Mesh::mesh msh;

    /** electrostatic sub-problem */
    electrostatSolver elec;

    /** this vector contains the material parameters for all regions for all the tetrahedrons */
    std::vector<Tetra::prm> paramTetra;

    /** this vector contains the material parameters for all surface regions for all the triangular facettes */
    std::vector<Facette::prm> paramFacette;

    /** if verbose set to true, some printing are sent to terminal */
    const bool verbose;

    /** maximum number of iteration for biconjugate stabilized gradient */
    const unsigned int MAXITER;

    /** number of Nodes (needed for templates) */
    const int NOD;

    /** returns Js = Ms/nu_0 */
    double getJs(Tetra::Tet const &tet) const;

    /** returns sigma of the tetraedron, (conductivity in (Ohm.m)^-1 */
    double getSigma(Tetra::Tet const &tet) const;
    
    /** \f$ \beta \f$ is polarization rate of the current */
    double getBeta(Tetra::Tet &tet) const;

    /** density of states at Fermi level, units : J^-1 nm^-3  */
    double getN0(Tetra::Tet &tet) const;

    /** length s-d */
    double getLsd(Tetra::Tet &tet) const;

    /** spin flip length */
    double getLsf(Tetra::Tet &tet) const;

    /** getter */
    Eigen::Vector3d get_Qn(Facette::Fac const &fac) const;
    
    /** affect extraField function and extraCoeffs_BE function for all the tetrahedrons */
    void prepareExtras(std::vector<Tetra::Tet> &v_tet, electrostatSolver &elec);

    /** solver, using biconjugate stabilized gradient, with diagonal preconditionner and Dirichlet
     * boundary conditions */
    int solve(const double _tol /**< [in] tolerance */ );

    std::vector<Eigen::Vector3d> Qs;

    static const int DIM_PROBLEM = 3;

    template <int N>
    void assembling(std::vector<int> &ind,//int ind[NBN],
                    Eigen::Matrix<double,DIM_PROBLEM*N,DIM_PROBLEM*N> &Ke,
                    Eigen::Matrix<double,DIM_PROBLEM*N,1> &Le,
                    algebra::w_sparseMat &K, std::vector <double> &L)
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
            for (int di=0; di<DIM_PROBLEM; di++) { L[di*NOD+i] += Le[di*N+ie]; }
            }
        }
    
    void integrales(Tetra::Tet &tet,
                    Eigen::Matrix<double,DIM_PROBLEM*Tetra::N,DIM_PROBLEM*Tetra::N> &AE,
                    Eigen::Matrix<double,DIM_PROBLEM*Tetra::N,1> &BE);
    
    void integrales(Facette::Fac &fac,
                    Eigen::Matrix<double,DIM_PROBLEM*Facette::N,1> &BE);
    };

#endif
