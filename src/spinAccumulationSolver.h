#ifndef spinAccumulationSolver_h
#define spinAccumulationSolver_h

#include <memory>
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
 slonczewskiInterface.h. If a tunnel barrier (tunnel) is given instead of jn, the normal current
 density depends on the magnetization (Robin condition for the potential, solved again at each
 call of compute()) and the spin is injected as a flux, see tunnelJunction.h. These
 surfaces are the ones given to the potential solver for its own boundary conditions. A surface may
 have both jn and s (s is then imposed on it).

 The spin accumulation s is continuous between the volume regions, as is the spin current. With
 s = mu_B N0 mu_s/2 (mu_s = (mu+ - mu-)/2, validation/valet_fert_1d), this is the continuity of the
 spin chemical potential mu_s of Valet-Fert if the regions have the same density of states N0
 (dens_state); otherwise, it imposes N0_A mu_s,A = N0_B mu_s,B on the interface between the regions
 A and B (validation/two_regions_1d). In verbose mode, a warning lists the neighbouring regions of
 different dens_state (checkDensityOfStates()).
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

    /** potential voltage(t) of F1 [V] at the time of the solve, for the tunnel barriers with
     * `V: voltage` (set by compute_all) */
    double voltage = NAN;

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

    /** anisotropic spin diffusion and longitudinal relaxation in the magnetic regions (settings
     * spin_accumulation.anisotropic_diffusion, default true), see integrales() */
    const bool anisotropicDiffusion;

    /** in verbose mode, warns about the neighbouring volume regions (sharing a face) of different
     * dens_state: s being continuous, mu_s is not continuous between them */
    void checkDensityOfStates(void) const;

    /** converts the polarizer stacks of the settings (physical units) with the density of states
     * of the volume region bordering each surface; checks that a surface with a N/F2 interface is
     * magnetic */
    void makeStacks(void);

    /** Robin boundary condition of the triangle f, on a surface with a polarizer stack
     * (slonczewski): adds the term in s to the matrix K and the source to the RHS L_rhs */
    void stackBoundaryCondition(Triangle::Tri &f);

    /** flux boundary condition of the triangle f, on a surface with a tunnel barrier: adds the
     * spin current (mu_B/e) jn P_T(m0) to the RHS L_rhs, jn = G_T(m0) (V - V1) */
    void tunnelBoundaryCondition(Triangle::Tri &f);

    /** solver of the electrostatic potential, kept to solve it again at each call of compute()
     * if a surface has a tunnel barrier (null if spin accumulation is disabled) */
    std::unique_ptr<electrostatSolver> pot_solver;

    /** true if a surface has a tunnel barrier: V depends on the magnetization */
    bool tunnel = false;

    /** computes all contributions to matrix AE from tetrahedron tet (LHS)
     * all = non magnetic metal + magnetic metal */
    void integrales(const Tetra::Tet &tet, Eigen::Matrix<double,DIM_PB*Tetra::N,DIM_PB*Tetra::N> &AE) const;

    /** computes magnetic metal contributions to spin diffusion from tetrahedron tet (RHS) */
    void integrales(const Tetra::Tet &tet /**< [in] */,
                    std::vector<double> &BE /**< [out] */) const;
    };

#endif
