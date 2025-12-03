#ifndef spinAccumulationSolver_h
#define spinAccumulationSolver_h

#include <vector>
#include "config.h"
#include "node.h"
#include "tetra.h"
#include "electrostatSolver.h"
#include "solver.h"
#include "meshUtils.h"

/** dimensionnality of the spin diffusion problem */
const int DIM_PB_SPIN_ACC = 3;

/** \class spinAcc
 container for Spin Accumulation constants, solver and related datas. The model obeys the diffusive
 spin equation. The boundary conditions are Dirichlet type. User has to provide through mesh and
 settings a surface S0 where spin diffusion vector s is a constant (zero recommended), and another
 surface S1 where current density and spin polarization is defined. This surface S1 must be the same
 as the one given to the potential solver for its own boundary conditions.
 */
class spinAcc : public solver<DIM_PB_SPIN_ACC>
    {
    public:
    /** constructor */
    spinAcc(Mesh::mesh &_msh /**< [in] ref to the mesh */,
    std::vector<Tetra::prm> & _pTetra /**< [in] ref to vector of param tetra (volume region parameters) */,
    std::vector<Facette::prm> & _pFac /**< [in] ref to vector of param facette (surface region parameters) */,
    const double _tol /**< [in] tolerance for bicg_dir solver */,  // _tol could be 1e-6
    const bool v /**< [in] verbose bool */,
    const int max_iter /**< [in] maximum number of iterations */):
        solver<DIM_PB_SPIN_ACC>(_msh,_pTetra,_pFac,"bicg_dir",_tol,v,max_iter), verbose(v), NOD(_msh.getNbNodes())
        {
        valDirichlet.resize(DIM_PROBLEM*NOD);
        boundaryConditions();
        }

    /** boundary conditions: a surface with a fixed s, and another surface with fixed normal current
     * density J and polarization vector P */
    void boundaryConditions(void);

    /** initializations: compute gradV and Hm and call prepareExtras method */
    void preCompute(void);

    /** call solver and update spin diffusion solution, returns true if solver succeeded */
    bool compute(void);

    /** set potential V */
    void setPotential(std::vector<double> &_V);

    /** check boundary conditions: mesh and settings have to define a single surface with constant
     * normal current density J, a vector polarization P and another single surface where spin diffusion = 0 */
    bool checkBoundaryConditions(bool verbose) const;

    /** spin accumulation solution over the nodes */
    std::vector<Eigen::Vector3d> s;

    private:
    /** container for potential values */
    std::vector<double> V;

    /** Dirichlet values of the components of s on the nodes, it is zero if the node is not in idxDirichlet */
    std::vector<double> valDirichlet;

    /** list of the indices for Dirichlet boundary conditions, it contains the indices of the nodes
     * where s is given by the user to be a constant on a surface */
    std::vector<int> idxDirichlet;

    /** fill valDirichlet and idxDirichlet vectors with k, a node index from fac.ind and the
     * corresponding spin diffusion s value*/
    void fillDirichletData(const int k, Eigen::Vector3d &s_value);

    /** table of the gradients of the potential, gradV.size() is the number of tetrahedrons
     * [gradV] = Volt/m = kg m A^-1 s^-3
     * */
    std::vector< Eigen::Matrix<double,Nodes::DIM,Tetra::NPI> > gradV;

    /** table of the Hm vectors (contribution of spinAcc to the tet::integrales) ; Hm.size() is the
     * number of tetra */
    std::vector< Eigen::Matrix<double,Nodes::DIM,Tetra::NPI> > Hm;

    /** if verbose set to true, some printing are sent to terminal */
    const bool verbose;

    /** number of digits in the optional output file */
    const int precision = 8;

    /** number of Nodes (needed for templates) */
    const int NOD;

    /** returns Js = mu_0 Ms */
    double getJs(Tetra::Tet const &tet) const;

    /** returns sigma of the tetraedron, (conductivity in (Ohm.m)^-1 */
    double getSigma(Tetra::Tet const &tet) const;
    
    /** \f$ P \f$ is polarization rate of the current */
    double getPolarization(Tetra::Tet &tet) const;

    /** density of states at Fermi level, units : J^-1 m^-3  */
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

    /** computes magnetic contributions to spin diffusion from tetrahedron tet */
    void integraleMag(Tetra::Tet &tet,
                      Eigen::Matrix<double,DIM_PROBLEM*Tetra::N,DIM_PROBLEM*Tetra::N> &AE,
                      std::vector<double> &BE);

    /** computes spin Hall effect contribution to spin diffusion. This contribution does not depend
     * on magnetization, its origin is Spin Orbit Torque. */
    void integraleSpinHall(Tetra::Tet &tet, std::vector<double> &BE);

    /** computes normal metal contributions to matrix AE and vector BE from tetrahedron tet */
    void integrales(Tetra::Tet &tet,
                    Eigen::Matrix<double,DIM_PROBLEM*Tetra::N,DIM_PROBLEM*Tetra::N> &AE);

    /** computes contributions to vector BE from facette fac */
    void integrales(Facette::Fac &fac, Eigen::Vector3d &Q, std::vector<double> &BE);
    };

#endif
