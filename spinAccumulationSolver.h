#ifndef spinAccumulationSolver_h
#define spinAccumulationSolver_h

#include <vector>
#include "config.h"
#include "node.h"
#include "tetra.h"
#include "electrostatSolver.h"
#include "solverUtils.h"

/** \class spinAcc
 container for Spin Accumulation constants, diffusive accumulation spin model, Boundary conditions
 (potential fixed value on one surface, current density on another surface)
 */

class spinAcc
    {
    public:
    /** constructor */
    spinAcc(Mesh::mesh &_msh /**< [in] ref to the mesh */,
    std::vector<Tetra::prm> _pTetra /**< [in] ref to vector of param tetra (volume region parameters) */,
    std::vector<Facette::prm> _pFac /**< [in] ref to vector of param facette (surface region parameters) */,
    const double _tol /**< [in] tolerance for solvers */,  // _tol could be 1e-6
    const bool v /**< [in] verbose bool */,
    const int max_iter /**< [in] maximum number of iteration */):
        msh(_msh), paramTetra(_pTetra),  paramFacette(_pFac),
        iter("bicg_dir",_tol,v,max_iter), verbose(v), NOD(_msh.getNbNodes()) {}

    /** boundary conditions on different surfaces. They must not share any triangle nor nodes. */
    Mesh::allBoundCond<Eigen::Vector3d> all_bc;

    /** initializations: compute gradV and Hm and call prepareExtras method */
    void preCompute(void);

    /** call solver and update spin diffusion solution, returns true if solver succeeded */
    bool compute(void);

    /** set potential V */
    void setPotential(std::vector<double> &_V);

    private:
    /** mesh object to store nodes, fac, tet, and others geometrical values related to the mesh */
    Mesh::mesh msh;

    /** this vector contains the material parameters for all regions for all the tetrahedrons */
    std::vector<Tetra::prm> paramTetra;

    /** this vector contains the material parameters for all surface regions for all the triangular facettes */
    std::vector<Facette::prm> paramFacette;

    /** container for potential values */
    std::vector<double> V;

    /** table of the gradients of the potential, gradV.size() is the number of tetra */
    std::vector< Eigen::Matrix<double,Nodes::DIM,Tetra::NPI> > gradV;

    /** table of the Hm vectors (contribution of spinAcc to the tet::integrales) ; Hm.size() is the
     * number of tetra */
    std::vector< Eigen::Matrix<double,Nodes::DIM,Tetra::NPI> > Hm;

    /** monitor the solver called in method solve() */
    algebra::iteration<double> iter;

    /** if verbose set to true, some printing are sent to terminal */
    const bool verbose;

    /** number of digits in the optional output file */
    const int precision = 8;

    /** number of Nodes (needed for templates) */
    const int NOD;

    /** returns Js = Ms/nu_0 */
    double getJs(Tetra::Tet const &tet) const;

    /** returns sigma of the tetraedron, (conductivity in (Ohm.m)^-1 */
    double getSigma(Tetra::Tet const &tet) const;
    
    /** \f$ P \f$ is polarization rate of the current */
    double getPolarization(Tetra::Tet &tet) const;

    /** density of states at Fermi level, units : J^-1 nm^-3  */
    double getN0(Tetra::Tet &tet) const;

    /** length s-d : only in magnetic material */
    double getLsd(Tetra::Tet &tet) const;

    /** spin flip length : exists in both non magnetic and magnetic metals */
    double getLsf(Tetra::Tet &tet) const;

    /** spin Hall constant */
    double getSpinHall(Tetra::Tet &tet) const;

    /** affect extraField function and extraCoeffs_BE function using lambdas for all the tetrahedrons (functor) */
    void prepareExtras(void);

    /** computes the gradient(V) for tetra tet */
    void calc_gradV(Tetra::Tet const &tet, Eigen::Ref<Eigen::Matrix<double,Nodes::DIM,Tetra::NPI>> _gradV);

    /** computes Hm contributions for each npi for tetrahedron tet */
    void calc_Hm(Tetra::Tet const &tet, Eigen::Ref<Eigen::Matrix<double,Nodes::DIM,Tetra::NPI>> _gradV,
                 Eigen::Ref<Eigen::Matrix<double,Nodes::DIM,Tetra::NPI>> _Hm);

    /** solver, using biconjugate stabilized gradient, with diagonal preconditionner and Dirichlet
     * boundary conditions */
    bool solve(void);

    /** dimensionnality of the spin accumulation problem */
    static const int DIM_PROBLEM = 3;

    /** computes contributions to matrix AE and vector BE from tetrahedron tet */
    void integrales(Tetra::Tet &tet,
                    Eigen::Matrix<double,DIM_PROBLEM*Tetra::N,DIM_PROBLEM*Tetra::N> &AE,
                    std::vector<double> &BE);

    /** computes contributions to vector BE from facette fac */
    void integrales(Facette::Fac &fac, Eigen::Vector3d &Q, std::vector<double> &BE);
    };

#endif
