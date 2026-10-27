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
 spin equation. User has to provide through mesh and settings at least one surface where spin
 diffusion vector s is a constant (zero recommended): Dirichlet boundary condition, and at least
 one surface where the current density jn is defined; if a spin polarization uP is also given on
 such a surface (and no s), the spin injection is a flux (Neumann) boundary condition. If a polarizer
 stack (slonczewski) is given instead, the injection is a Robin boundary condition, see
 slonczewskiInterface.h. These
 surfaces are the ones given to the potential solver for its own boundary conditions. A surface may
 have both jn and s (s is then imposed on it).
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
            const int max_iter /**< [in] maximum number of iterations */);

    /** Dirichlet boundary conditions: all the surfaces with a fixed s (the surface with fixed normal
     * current density J and polarization vector P, without s, is a flux condition, added to the RHS
     * in solve())
     * set valDirichlet values and fill vector of indices idxDirichlet
     * */
    void boundaryConditions(void); // should be private

    /** call solver and update spin diffusion solution, returns true if solver succeeded */
    bool compute(void);

    /** dimensionless factor of the electric sources at the time of the solve:
     * current_density(t)/j_ref if a surface region has `jn: current_density`, its jn being the
     * reference j_ref = 1 A/m^2, 1 otherwise. The electrokinetics being linear, the potential V0
     * is solved once with the jn of the settings, and grad V = jn_factor grad V0,
     * jn -> jn_factor jn (the constant imposed V does not change grad V). It multiplies the
     * source terms (polarized current, injection through the surface with jn and uP), not the
     * imposed s */
    double jn_factor = 1.0;

    /** electrostatic potential on the nodes at the time of the last solve [V], indexed as s:
     * V_D + jn_factor V, V being solved with V = 0 on the surface with V and V_D the imposed V
     * (the electrokinetics being linear, the potential at time t follows the current density;
     * jn_factor = 1 for a numeric jn) */
    std::vector<double> potential(void) const;

    /** solution of the spin diffusion vector s over the nodes */
    std::vector<Eigen::Vector3d> s;

    /** check boundary conditions: mesh and settings have to define at least one surface with a
     * normal current density jn and at least one surface where the spin diffusion s is given */
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

    /** electrostatic potential on the nodes, solved once with V = 0 on the surface with V (V_D
     * is added by potential()), for the jn of the settings, or, with `jn: current_density`, for
     * the reference j_ref = 1 A/m^2 (V_ref); grad V = jn_factor grad V at time t, jn_factor being
     * dimensionless
     * unit [V] = Volt = kg A^-1 s^-3 */
    std::vector<double> V;

    /** potential imposed on the surface with V (Dirichlet condition of the electrostatic
     * problem) [V] */
    double V_D = 0.0;

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

    /** nodal magnetization of the tetrahedron (NEXT step, the latest computed one) used by the
     * spin diffusion problem, columns are the nodes */
    Eigen::Matrix<double,Nodes::DIM,Tetra::N> calc_u_nod(const Tetra::Tet &tet) const;

    /** affect extraField member function of all tetrahedrons
     * extraField is computing the contribution from the spin diffusion s to llg
     * */
    void prepareExtraField(void) const;

    /** solver, using biconjugate stabilized gradient, with diagonal preconditionner and Dirichlet
     * boundary conditions */
    bool solve(void);

    /** polarizer stacks in the units of the spin diffusion problem, indexed as paramTri (only the
     * surfaces with slonczewski are set) */
    std::vector<Slonczewski::Stack> stacks;

    /** converts the polarizer stacks of the settings (physical units) with the density of states
     * of the volume region bordering each surface; checks that a surface with a N/F2 interface is
     * magnetic */
    void makeStacks(void);

    /** Robin boundary condition of the triangle f, on a surface with a polarizer stack
     * (slonczewski): adds the term in s to the matrix K and the source to the RHS L_rhs */
    void stackBoundaryCondition(Triangle::Tri &f);

    /** computes all contributions to matrix AE from tetrahedron tet (LHS)
     * all = non magnetic metal + magnetic metal */
    void integrales(const Tetra::Tet &tet, Eigen::Matrix<double,DIM_PB*Tetra::N,DIM_PB*Tetra::N> &AE) const;

    /** computes magnetic metal contributions to spin diffusion from tetrahedron tet (RHS) */
    void integrales(const Tetra::Tet &tet /**< [in] */,
                    std::vector<double> &BE /**< [out] */) const;
    };

#endif
