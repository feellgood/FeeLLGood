/** \file electrostatSolver.h
  \brief solver for electrostatic problem when spin accumulation is required
  header containing electrostatSolver class. It uses biconjugate stabilized gradient with diagonal
  preconditioner. The solver is only called once to compute voltages V for each nodes of the mesh,
  when spin accumulation computation is involved.
 */

#include <iostream>
#include <map>

#include "config.h"
#include "fem.h"
#include "mesh.h"
#include "algebra/algebra.h"

/** \class electrostatSolver
this class is containing both data and a solver to compute potential from dirichlet boundary
conditions problem for the current density flowing in the sample.
*/
class electrostatSolver
    {
public:
    /** constructor */
    electrostatSolver(
            Mesh::mesh & _msh /**< [in] reference to the mesh */,
            std::vector<Tetra::prm> _pTetra /**< [in] ref to vector of param tetra (volume region parameters) */,
            std::vector<Facette::prm> _pFac /**< [in] ref to vector of param facette (surface region parameters) */,
            const double _tol /**< [in] tolerance for solvers */,
            const bool v /**< [in] verbose mode for iteration monitor */,
            const int max_iter /**< [in] maximum number of iteration */)
        : msh(_msh), paramTetra(_pTetra), paramFacette(_pFac), iter("cg_dir",_tol,v,max_iter) {}

    /** electrostatic potential values for boundary conditions, V.size() is the size of the vector
     * of nodes */
    std::vector<double> V;

    /** check boundary conditions: mesh and settings have to define a single surface with constant current
     * density J (normal componant to the surface) and another single surface with constant potential V */
    bool checkBoundaryConditions(void) const;

    /** basic informations on boundary conditions */
    void infos(void);

    /** compute side problem (electrostatic potential on the nodes) integrales for matrix
     * coefficients,inputs from tet */
    void integrales(Tetra::Tet const &tet, Eigen::Ref<Eigen::Matrix<double,Tetra::N,Tetra::N> > AE);

    /** compute integrales for vector coefficients, input from facette */
    void integrales(Facette::Fac const &fac, std::vector<double> &BE);

    /** text file (tsv) writing function for the solution V over all volume regions of the mesh,
     * node indices are zero based */
    bool save(const std::string V_fileName /**< [in] output file name */, std::string const &metadata /**< [in] */) const;

    /** returns sigma of the tetraedron, (conductivity in (Ohm.m)^-1 */
    double getSigma(Tetra::Tet const &tet) const;

    /** returns current density of the facette if it is defined in the boundary conditions, else zero */
    double getCurrentDensity(Facette::Fac const &fac) const;

    /** solves the potential and stores result in V, save to text file if needed
     * if verbose set to true, some printing are sent to terminal */
    void compute(const bool verbose /**< [in] */, const std::string V_fileName /**< [in] output file name for V solution */);
private:
    /** mesh object to store nodes, fac, tet, and others geometrical values related to the mesh (
     * const ref ) */
    Mesh::mesh msh;

    /** this vector contains the material parameters for all volume regions for all the tetrahedrons */
    std::vector<Tetra::prm> paramTetra;

    /** this vector contains the material parameters for all surface regions for all the triangular facettes */
    std::vector<Facette::prm> paramFacette;

    /** monitor the solver called in method solve() */
    algebra::iteration<double> iter;

    /** number of digits in the optional output file */
    const int precision = 8;

    /** solver, using conjugate gradient with masking, with diagonal preconditionner and Dirichlet
     * boundary conditions, returns true if has converged */
    bool solve(void);
    };  // end class electrostatSolver
