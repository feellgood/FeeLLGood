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
            const double _sigma /**< [in] conductivity */,
            const double _tol /**< [in] tolerance for solvers */,
            const bool v /**< [in] verbose bool */,
            const int max_iter /**< [in] maximum number of iteration */,
            const std::string _fileName /**< [in] output .sol file name for electrostatic potential */)
        : msh(_msh), sigma(_sigma), verbose(v), MAXITER(max_iter), fileName(_fileName)
        {
        if (verbose)
            { infos(); }
        V.resize(_msh.getNbNodes());
        bool has_converged = solve(_tol);
        if (has_converged)
            {
            if (V_file)
                {
                if (verbose)
                    { std::cout << "writing electrostatic potential to file " << fileName << std::endl; }
                bool iznogood = save("## columns: index\tV\n");
                if (verbose && iznogood)
                    { std::cout << "file " << fileName << " written.\n"; }
                }
            std::for_each(msh.tet.begin(), msh.tet.end(), [this](Tetra::Tet const &tet)
                          {
                          Eigen::Matrix<double,Nodes::DIM,Tetra::NPI> _gradV;
                          calc_gradV(tet, _gradV);
                          gradV.push_back(_gradV);

                          Eigen::Matrix<double,Nodes::DIM,Tetra::NPI> _Hm;
                          calc_Hm(tet, _gradV, _Hm);
                          Hm.push_back(_Hm);
                          });
            }
        else
            {
            std::cerr << "Solver (ElectroStatic) has not converged" << std::endl;
            exit(1);
            }
        }

    /** computes the gradient(V) for tetra tet */
    inline void calc_gradV(Tetra::Tet const &tet, Eigen::Ref<Eigen::Matrix<double,Nodes::DIM,Tetra::NPI>> _gradV)
        {
        for (int npi = 0; npi < Tetra::NPI; npi++)
            {
            Eigen::Vector3d v(0,0,0);
            for (int i = 0; i < Tetra::N; i++)
                {
                v += V[tet.ind[i]] * tet.da.row(i);
                }
            _gradV.col(npi) = v;
            }
        }

    /** computes Hm contributions for each npi for tetrahedron tet */
    inline void calc_Hm(Tetra::Tet const &tet, Eigen::Ref<Eigen::Matrix<double,Nodes::DIM,Tetra::NPI>> _gradV,
                 Eigen::Ref<Eigen::Matrix<double,Nodes::DIM,Tetra::NPI>> _Hm)
        {
        Eigen::Matrix<double,Nodes::DIM,Tetra::NPI> p_g;
        tet.getPtGauss(p_g);

        for (int npi = 0; npi < Tetra::NPI; npi++)
            { _Hm.col(npi) = -sigma * _gradV.col(npi).cross(p_g.col(npi)); }
        }



    /** electrostatic potential values for boundary conditions, V.size() is the size of the vector
     * of nodes */
    std::vector<double> V;

    /** table of the gradients of the potential, gradV.size() is the number of tetra */
    std::vector< Eigen::Matrix<double,Nodes::DIM,Tetra::NPI> > gradV;

    /** table of the Hm vectors (contribution of spinAcc to the tet::integrales) ; Hm.size() is the
     * number of tetra */
    std::vector< Eigen::Matrix<double,Nodes::DIM,Tetra::NPI> > Hm;

    /** basic informations on boundary conditions */
    void infos(void);

    /** assemble the matrix K from tet and Ke inputs */
    void assembling_mat(Tetra::Tet const &tet, double Ke[Tetra::N][Tetra::N], std::vector<Eigen::Triplet<double>> &K);

    /** assemble the vector L from fac and Le inputs */
    void assembling_vect(Facette::Fac const &fac, std::vector<double> const &Le, Eigen::Ref<Eigen::VectorXd> L);

    /** compute side problem (electrostatic potential on the nodes) integrales for matrix
     * coefficients,input from tet ; sigma is the region conductivity */
    void integrales(Tetra::Tet const &tet, double sigma, double AE[Tetra::N][Tetra::N]);

    /** compute integrales for vector coefficients, input from facette */
    void integrales(Facette::Fac const &fac, double pot_val, std::vector<double> &BE);

    /** text file (tsv) writing function for the solution V over all volume regions of the mesh,
     * node indices are zero based */
    bool save(std::string const &metadata /**< [in] */) const;

private:
    /** mesh object to store nodes, fac, tet, and others geometrical values related to the mesh (
     * const ref ) */
    Mesh::mesh msh;

    /** Conductivity Ohm^-1 nm^-1 */
    double sigma;

    /** if verbose set to true, some printing are sent to terminal */
    const bool verbose;

    /** maximum number of iteration for biconjugate stabilized gradient */
    const unsigned int MAXITER;  // fixed to 5000 in ref code

    /** number of digits in the optional output file */
    const int precision = 8;

    /** if true a text file V.sol is generated. It contains the solution of the
                     electrostatic problem on the nodes of the mesh */
    const bool V_file = false;

    /** output .sol file name for electrostatic problem */
    const std::string fileName;

    /** boundary conditions, stored as a vector of pairs.
    First element of the pair is the surface region name given in the mesh by its physical name;
    Second element is the electrostatic potential associated value  */
    std::vector<std::pair<std::string, double>> boundaryCond;

    /** fill matrix and vector to solve potential values on each node */
    void prepareData(std::vector<Eigen::Triplet<double>> &Kw, Eigen::Ref<Eigen::VectorXd> Lw);

    /** solver, using biconjugate stabilized gradient, with diagonal preconditionner and Dirichlet
     * boundary conditions */
    int solve(const double _tol);

    };  // end class electrostatSolver
