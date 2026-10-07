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
 diffusion vector s is a constant (zero recommended): Dirichlet boundary condition, and one surface
 S1 where the current density jn is defined; if a spin polarization uP is also given on S1 (and no
 s), the spin injection is a flux (Neumann) boundary condition. This surface S1 must be the same
 as the one given to the potential solver for its own boundary conditions. A surface may have both
 jn and s (s is then imposed on it).
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

    /** solution of the spin diffusion vector s over the nodes */
    std::vector<Eigen::Vector3d> s;

    /** check boundary conditions: mesh and settings have to define a single surface with constant
     * normal current density J and at least one surface where the spin diffusion s is given */
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

    /** electrostatic potential V on the nodes
     * unit [V] = Volt = kg A^-1 s^-3 */
    std::vector<double> V;

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

    /** computes all contributions to matrix AE from tetrahedron tet (LHS)
     * all = non magnetic metal + magnetic metal */
    void integrales(const Tetra::Tet &tet, Eigen::Matrix<double,DIM_PB*Tetra::N,DIM_PB*Tetra::N> &AE) const;

    /** computes magnetic metal contributions to spin diffusion from tetrahedron tet (RHS) */
    void integrales(const Tetra::Tet &tet /**< [in] */,
                    std::vector<double> &BE /**< [out] */) const;
    };

#endif
