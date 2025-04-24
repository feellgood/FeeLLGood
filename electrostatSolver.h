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

/** dimensionnality of the problem to solve.
 * here electrostatic problem to solve computes scalar potential V on the nodes */
const int DIM_PROBLEM = 1;

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
            const bool v /**< [in] verbose bool */,
            const int max_iter /**< [in] maximum number of iteration */,
            const bool _V_file /**< [in] if true an output tsv file containing potential V is written */,
            const std::string _fileName /**< [in] output .sol file name for electrostatic potential */)
        : msh(_msh), paramTetra(_pTetra), paramFacette(_pFac), verbose(v), MAXITER(max_iter), V_file(_V_file), fileName(_fileName), NOD(_msh.getNbNodes())
        {
        if (verbose)
            { infos(); }
        V.resize(NOD);
        bool has_converged = solve(_tol);
        if (has_converged)
            {
            std::cout << "Solver (ElectroStatic) has converged.\n";
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
            std::cerr << "Solver (ElectroStatic) has not converged.\n";
            exit(1);
            }
        }

    /** computes the gradient(V) for tetra tet */
    void calc_gradV(Tetra::Tet const &tet, Eigen::Ref<Eigen::Matrix<double,Nodes::DIM,Tetra::NPI>> _gradV);

    /** computes Hm contributions for each npi for tetrahedron tet */
    void calc_Hm(Tetra::Tet const &tet, Eigen::Ref<Eigen::Matrix<double,Nodes::DIM,Tetra::NPI>> _gradV,
                 Eigen::Ref<Eigen::Matrix<double,Nodes::DIM,Tetra::NPI>> _Hm);

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

    /** compute side problem (electrostatic potential on the nodes) integrales for matrix
     * coefficients,inputs from tet */
    void integrales(Tetra::Tet const &tet, Eigen::Ref<Eigen::Matrix<double,Tetra::N,Tetra::N> > AE);

    /** compute integrales for vector coefficients, input from facette */
    void integrales(Facette::Fac const &fac, Eigen::Ref<Eigen::Matrix<double,Facette::N,1> > BE);

    /** text file (tsv) writing function for the solution V over all volume regions of the mesh,
     * node indices are zero based */
    bool save(std::string const &metadata /**< [in] */) const;

    /** returns sigma of the tetraedron, (conductivity in (Ohm.m)^-1 */
    double getSigma(Tetra::Tet const &tet) const;

    /** returns current density of the facette if it is defined in the boundary conditions, else zero */
    double getCurrentDensity(Facette::Fac const &fac) const;

private:
    /** mesh object to store nodes, fac, tet, and others geometrical values related to the mesh (
     * const ref ) */
    Mesh::mesh msh;

    /** this vector contains the material parameters for all volume regions for all the tetrahedrons */
    std::vector<Tetra::prm> paramTetra;

    /** this vector contains the material parameters for all surface regions for all the triangular facettes */
    std::vector<Facette::prm> paramFacette;

    /** if verbose set to true, some printing are sent to terminal */
    const bool verbose;

    /** maximum number of iteration for biconjugate stabilized gradient */
    const unsigned int MAXITER;  // fixed to 5000 in ref code

    /** number of digits in the optional output file */
    const int precision = 8;

    /** if true a text file V.sol is generated. It contains the solution of the
                     electrostatic problem on the nodes of the mesh */
    const bool V_file;

    /** output .sol file name for electrostatic problem */
    const std::string fileName;

    /** number of Nodes (needed for templates) */
    const int NOD;

/** assemble the matrix K from Ke input */
template <int N>
    void assemble_mat(std::vector<int> &ind,
                      Eigen::Matrix<double,DIM_PROBLEM*N,DIM_PROBLEM*N> &Ke, algebra::w_sparseMat &K)
        {
        for (int ie=0; ie<N; ie++)
            {
            int i= ind[ie];
            for (int je=0; je<N; je++)
                {
                int j= ind[je];
                for (int di=0; di<DIM_PROBLEM; di++)
                    for (int dj=0; dj<DIM_PROBLEM; dj++)
                        K.insert(di*NOD+i, dj*NOD+j, Ke(di*N+ie,dj*N+je));
	            }
            }
        }

/** assemble the vector L from Le input */
template <int N>
    void assemble_vect(std::vector<int> &ind,
                       Eigen::Matrix<double,DIM_PROBLEM*N,1> &Le, std::vector <double> &L)
        {
        for (int ie=0; ie<N; ie++)
            {
            int i= ind[ie];
            for (int di=0; di<DIM_PROBLEM; di++) { L[di*NOD+i] += Le[di*N+ie]; }
            }
        }

/**
 suppress twins in the indices Dirichlet list v_idx
*/
template <typename T>
void suppress_copies(std::vector<T> &v_idx)
    {
    std::sort (v_idx.begin(), v_idx.begin()+v_idx.size());
    auto it = std::unique (v_idx.begin(), v_idx.end() );
    v_idx.resize( std::distance(v_idx.begin(),it) );
    v_idx.shrink_to_fit();
    }

    /** solver, using biconjugate stabilized gradient, with diagonal preconditionner and Dirichlet
     * boundary conditions */
    int solve(const double _tol);
    };  // end class electrostatSolver
