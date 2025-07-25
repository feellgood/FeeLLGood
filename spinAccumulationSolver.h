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

class spinAcc
    {
    public:
    /** constructor */
    spinAcc(Mesh::mesh &_msh /**< [in] ref to the mesh */,
    electrostatSolver &_elec /**< [in] ref to the the electrostatic sub_problem */,
    std::vector<Tetra::prm> _pTetra /**< [in] ref to vector of param tetra (volume region parameters) */,
    std::vector<Facette::prm> _pFac /**< [in] ref to vector of param facette (surface region parameters) */,
    const double _tol /**< [in] tolerance for solvers */,  // _tol could be 1e-6
    const bool v /**< [in] verbose bool */,
    const int max_iter /**< [in] maximum number of iteration */):
        msh(_msh), elec(_elec), paramTetra(_pTetra),  paramFacette(_pFac),
        iter("bicg_dir",_tol,v,max_iter), verbose(v), NOD(_msh.getNbNodes())
        {
        Qs.resize(NOD);
        bool has_converged = solve();
        if (!has_converged)
            { std::cout << "spin accumulation solver: " << iter.infos() << std::endl; exit(1); }
        prepareExtras();
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

    /** monitor the solver called in method solve() */
    algebra::iteration<double> iter;

    /** if verbose set to true, some printing are sent to terminal */
    const bool verbose;

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
    
    /** affect extraField function and extraCoeffs_BE function using lambdas for all the tetrahedrons (functor) */
    void prepareExtras(void);

    /** solver, using biconjugate stabilized gradient, with diagonal preconditionner and Dirichlet
     * boundary conditions */
    int solve(void);

    /** solution of the accumulation spin diffusion problem (3D vector field) */
    std::vector<Eigen::Vector3d> Qs;

    /** dimensionnality of the spin accumulation problem */
    static const int DIM_PROBLEM = 3;

    /** computes contributions to matrix AE and vector BE from tetrahedron tet */
    void integrales(Tetra::Tet &tet,
                    Eigen::Matrix<double,DIM_PROBLEM*Tetra::N,DIM_PROBLEM*Tetra::N> &AE,
                    std::vector<double> &BE);

    /** computes contributions to vector BE from facette fac */
    void integrales(Facette::Fac &fac, std::vector<double> &BE);
    };

#endif
