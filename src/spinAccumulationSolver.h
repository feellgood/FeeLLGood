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
    /** spin accumulation constructor
     * if mySettings.spin_acc bool is false it is a do nothing constructor
     * */
    spinAcc(const Settings &mySettings /**< [in] */,
            Mesh::mesh &_msh /**< [in] ref to the mesh */,
            const double _tol /**< [in] tolerance for bicg_dir solver */,
            const int max_iter /**< [in] maximum number of iterations */):
            solver<DIM_PB_SPIN_ACC>(_msh, mySettings.paramTetra, mySettings.paramTriangle,
                                    "bicg_dir", _tol, mySettings.verbose, max_iter)
        {
        if (mySettings.spin_acc)
            {
            checkBoundaryConditions();
            valDirichlet.resize(DIM_PB*NOD);
            boundaryConditions();
            electrostatSolver pot_solver(_msh, mySettings.paramTetra, mySettings.paramTriangle,
                                         1e-8, mySettings.verbose, 1000);
            pot_solver.checkBoundaryConditions();
            pot_solver.V.resize(_msh.getNbNodes());
            std::string V_fileName("");
            if(mySettings.V_file)
                V_fileName = mySettings.getSimName() + "_V.sol";
            pot_solver.compute(mySettings.verbose, V_fileName);
            preCompute(pot_solver.V);
            if(!compute())
                {
                std::cout << "Error: spin diffusion solver(first try) failed.\n";
                exit(1);
                }
            }
        }

    /** boundary conditions: a surface with a fixed s, and another surface with fixed normal current
     * density J and polarization vector P
     * set valDirichlet values and fill vector of indices idxDirichlet
     * */
    void boundaryConditions(void); // should be private

    /** initializations: compute gradV and call prepareExtras method */
    void preCompute(std::vector<double> &V);

    /** call solver and update spin diffusion solution, returns true if solver succeeded */
    bool compute(void);

    /** solution of the spin diffusion vector s over the nodes */
    std::vector<Eigen::Vector3d> s;

    /** check boundary conditions: mesh and settings have to define a single surface with constant
     * normal current density J, a vector polarization P and another single surface where spin
     * diffusion = 0 */
    void checkBoundaryConditions(void) const override;

    private:
    /** Dirichlet values of the components of s on the nodes, it is zero if the node is not in
     * idxDirichlet */
    std::vector<double> valDirichlet;

    /** list of the indices for Dirichlet boundary conditions, it contains the indices of the nodes
     * where s is given by the user to be a constant on a surface */
    std::vector<int> idxDirichlet;

    /** fill valDirichlet and idxDirichlet vectors with k, a node index from tri.ind and the
     * corresponding spin diffusion s value*/
    void fillDirichletData(const int k, Eigen::Vector3d &s_value);

    /** table of the gradients of the potential, gradV.size() is the number of tetrahedrons
     * [gradV] = Volt/m = kg m A^-1 s^-3
     * */
    std::vector< Eigen::Matrix<double,Nodes::DIM,Tetra::NPI> > gradV;

    /** number of digits in the optional output file */
    const int precision = 8;

    /** returns Ms */
    double getMs(const Tetra::Tet &tet) const;

    /** returns sigma of the tetraedron, (conductivity in (Ohm.m)^-1 */
    double getSigma(const Tetra::Tet &tet) const;

    /** \f$ P \f$ is polarization rate of the current density */
    double getPolarizationRate(const Tetra::Tet &tet) const;

    /** diffusion constant, units: s^-1 m^2 */
    double getDiffusionCst(const Tetra::Tet &tet) const;

    /** length s-d : only in magnetic material */
    double getLsd(const Tetra::Tet &tet) const;

    /** spin flip length : exists in both non-magnetic and magnetic metals */
    double getLsf(const Tetra::Tet &tet) const;

    /** affect extraField member function of all tetrahedrons */
    void prepareExtras(void);

    /** solver, using biconjugate stabilized gradient, with diagonal preconditionner and Dirichlet
     * boundary conditions */
    bool solve(void);

    /** computes all contributions to matrix AE from tetrahedron tet (LHS)
     * all = non magnetic metal + magnetic metal */
    void integrales(const Tetra::Tet &tet, Eigen::Matrix<double,DIM_PB*Tetra::N,DIM_PB*Tetra::N> &AE) const;

    /** computes magnetic metal contributions to spin diffusion from tetrahedron tet (RHS) */
    void integrales(Tetra::Tet &tet,std::vector<double> &BE);
    };

#endif
