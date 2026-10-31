#ifndef electrostat_h
#define electrostat_h

/** \file electrostatSolver.h
  \brief solver for electrostatic problem when spin accumulation is required
  header containing electrostatSolver class. It uses conjugate gradient with diagonal
  preconditioner. The solver is called once to compute voltages V for each nodes of the mesh,
  when spin accumulation computation is involved, or at each time step if a surface has a tunnel
  barrier, its conductance depending on the magnetization (tunnelJunction.h).
 */

#include "solver.h"
#include "config.h"
#include "algebra/algebra.h"

/** dimensionnality of the electrostatic problem */
constexpr int DIM_PB_ELEC = 1;

/** \class electrostatSolver
this class is containing both data and a solver to compute potential from dirichlet boundary
conditions problem for the current density flowing in the sample.
*/
class electrostatSolver: public solver<DIM_PB_ELEC>
    {
public:
    /** constructor */
    electrostatSolver(
            Mesh::mesh & _msh /**< [in] reference to the mesh */,
            const std::vector<Tetra::prm> & _pTetra /**< [in] ref to vector of param tetra
                                                (volume region parameters) */,
            const std::vector<Triangle::prm> & _pTri /**< [in] ref to vector of param triangle
                                                (surface region parameters) */,
            const double _tol /**< [in] cg_dir tolerance */,
            const bool v /**< [in] verbose mode for iteration monitor */,
            const int max_iter /**< [in] maximum number of iterations */):
            solver<DIM_PB_ELEC>(_msh, _pTetra, _pTri, "cg_dir", _tol, v, max_iter) {}

    /** electrostatic potential on the nodes, solved with V = 0 on the surface with V (the imposed
     * V is added afterwards), V.size() is the size of the vector of nodes */
    std::vector<double> V;

    /** potential voltage(t) of F1 [V] at the time of the solve, for the tunnel barriers with
     * `V: voltage` (set by spinAcc) */
    double voltage = NAN;

    /** basic informations on boundary conditions */
    void infos(void) const;

    /** compute side problem (electrostatic potential on the nodes) integrales for matrix
     * coefficients,inputs from tet */
    void integrales(const Tetra::Tet &tet, Eigen::Ref<Eigen::Matrix<double,Tetra::N,Tetra::N> > AE) const;

    /** compute integrales for vector coefficients, input from triangle */
    void integrales(const Triangle::Tri &tri, std::vector<double> &BE) const;

    /** compute the integrales of a triangle of a surface with a tunnel barrier, Robin condition
     * jn = G_T(m0) (V - V1): G_T a_i a_j to the matrix, G_T (V1 - V_D) a_i to the vector (V being
     * solved with V = 0 on the surface with V), m0 being the magnetization at the Gauss points */
    void integrales(const Triangle::Tri &tri,
                    Eigen::Ref<Eigen::Matrix<double,Triangle::N,Triangle::N> > AE,
                    std::vector<double> &BE) const;

    /** true if a surface region has a tunnel barrier: the potential depends on the magnetization
     * and has to be solved again when it changes */
    bool hasTunnel(void) const;

    /** text file (tsv) writing function for the solution V over all volume regions of the mesh,
     * node indices are zero based */
    bool save(const std::string &V_fileName /**< [in] output file name */,
              const std::string &metadata /**< [in] */) const;

    /** returns sigma of the tetraedron, (conductivity in (Ohm.m)^-1 */
    double getSigma(const Tetra::Tet &tet) const;

    /** returns current density of the triangle if it is defined in the boundary conditions,
     * else zero */
    double getCurrentDensity(const Triangle::Tri &tri) const;

    /** solves the potential and stores result in V, save to text file if needed
     * if verbose set to true, some printing are sent to terminal */
    void compute(bool verbose /**< [in] */,
            const std::string& V_fileName /**< [in] output file name for V solution */);

    /** check boundary conditions: mesh and settings have to define at least one surface with a
     * normal current density jn, either all numeric or all `jn: current_density`, or with a tunnel
     * barrier (not mixed with `jn: current_density`), and a single surface with constant
     * potential V */
    void checkBoundaryConditions(void) const override;

private:
    /** number of digits in the optional output file */
    const int precision = 8;

    /** potential imposed on the surface with V [V], added to V in the output file */
    double V_D = 0.0;

    /** true if a surface has jn: current_density: V is then the reference potential V_ref (for
     * j_ref = 1 A/m^2), written as is in the output file */
    bool reference_problem = false;

    /** solver, using conjugate gradient with masking, with diagonal preconditionner
     * and Dirichlet boundary conditions, returns true if has converged; the system is assembled
     * again at each call, and the solve starts from the current V */
    bool solve(void);
    };  // end class electrostatSolver

#endif
