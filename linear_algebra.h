#ifndef linear_algebra_h
#define linear_algebra_h

/** \file linear_algebra.h
\brief secondary header, it grabs altogether the linear algebra by the solver to apply fem method
<br> It encapsulates the calls to eigen BiCGSTAB solver, the assemblage and projection of the matrix for all elements
<br> projection and matrix assembly is multithreaded for tetrahedron, monothread for facette
*/
#include <random>
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wmaybe-uninitialized"
#include <execution>
#pragma GCC diagnostic pop

#include "config.h"

#include "facette.h"
#include "feellgoodSettings.h"
#include "mesh.h"
#include "node.h"
#include "tetra.h"

#include "algebra/algebra.h"

/** \class LinAlgebra
convenient class to grab altogether some part of the calculations involved using algebra::bicg solver at each
timestep. The solver is handled by solver method, and is using algebra sparse matrices(Row major). The write sparse matrix
is prepared in 'batch mode', to add all non zero coefficients with a '+=' logic. Then it is turned into a read sparse matrix before being used by bicg algo.
Be aware of time units: when entering solver method, division by gamma0 and multiplication by gamma0 when ending are mandatory.
The bicg algorithm is monitored by iter object.
When debugging it might be usefull to set iter verbosity differently from LinAlgebra Solver (see constructor list initialization)
*/
class LinAlgebra
    {
public:
    /** constructor
    When debugging it might be usefull to set iter verbosity differently
    */
    LinAlgebra(Settings &s /**< [in] */, Mesh::mesh &my_msh /**< [in] */)
        : NOD(my_msh.getNbNodes()), refMsh(&my_msh), iter("bicg",s.TOL,s.verbose,s.MAXITER),
          K(build_shape()), verbose(s.verbose),
          prmTetra(s.paramTetra), prmFacette(s.paramFacette)
        {
        L_rhs.resize(2*NOD);
        Xw.resize(2*NOD);
        base_projection();
        if (!s.recenter)
            { idx_dir = Nodes::IDX_UNDEF; }
        else
            { idx_dir = s.recentering_direction; }

        if(s.getFieldType() == R4toR3)
            { setExtSpaceField(s); }
        v_max = 0;
        my_msh.build_lvd(lvd);
        }

    /** build a matrix shape suitable for our mesh */
    algebra::MatrixShape build_shape()
        {
        algebra::MatrixShape shape(2 * NOD);
        auto insert_pair = [this, &shape](int i, int j)
            {
            shape[2*i].insert(2*j);
            shape[2*i].insert(2*j+1);
            shape[2*i+1].insert(2*j);
            shape[2*i+1].insert(2*j+1);
            };
        for (int i = 0; i < NOD; ++i)
            { insert_pair(i, i); }
        for (auto edge: refMsh->edges)
            {
            insert_pair(edge.first, edge.second);
            insert_pair(edge.second, edge.first);
            }
        return shape;
        }

    /** computes inner data structures of tetraedrons and triangular facettes (K matrices and L vectors)
    this member function is overloaded to fit to two different situations, either if std::function passed to element = tetra is corresponding to the simple case of constant external field applied to the magnetic region or space dependant. Here is the constant space applied field.
    */
    void prepareElements(Eigen::Vector3d const &Hext /**< [in] applied field */, timing const &t_prm /**< [in] */);

    /** computes inner data structures of tetraedrons and triangular facettes (K matrices and L vectors) 
    this member function is overloaded to fit to two different situations, either if std::function passed to element = tetra is corresponding to the simple case of constant external field applied to the magnetic region or space dependant. Here is the variable space applied field.
    */
    void prepareElements(double const A_Hext /**< [in] amplitude applied field (might be time dependant)*/, timing const &t_prm /**< [in] */);

    /** build init guess for bicg solver */
    void buildInitGuess(std::vector<double> &G/**< [out] */) const;

    /** solver, uses stabilized biconjugate gradient solver (bicg) with diagonal preconditionner, sparse matrix and vector are filled with multiThreading. Sparse matrix is row major.
    */
    int solver(timing const &t_prm /**< [in] */);

    /** setter for DW_dz */
    inline void set_DW_vz(double vz /**< [in] */) { DW_vz = vz; }

    /** getter for v_max */
    inline double get_v_max(void) { return v_max; }

    /** when external applied field is of field_type R4toR3 values of field_space are stored in spaceField */
    void setExtSpaceField(Settings &s /**< [in] */);

    /** computes local vector basis {ep,eq} in the tangeant plane for projection on the elements */
    void base_projection();
private:
    /** recentering index direction if any */
    Nodes::index idx_dir;

    /** number of nodes, also an offset for filling sparseMatrix, initialized by constructor */
    const int NOD;

    /** direct access to the mesh */
    Mesh::mesh *refMsh;

    /** monitor the solver */
    algebra::iteration<double> iter;

    /** matrix of the system to solve */
    algebra::r_sparseMat K;

    /** RHS vector of the system to solve */
    std::vector<double> L_rhs;

    /** solution of the system to solve */
    std::vector<double> Xw;

    /** verbosity */
    const int verbose;

    /** material parameters of the tetrahedrons */
    const std::vector<Tetra::prm> &prmTetra;

    /** material parameters of the facettes */
    const std::vector<Facette::prm> &prmFacette;

    /** speed of the domain wall */
    double DW_vz;

    /** maximum speed of the magnetization in the whole physical object */
    double v_max;

    /** external applied space field, values on gauss points, size is number of tetraedrons */
    std::vector< Eigen::Matrix<double,Nodes::DIM,Tetra::NPI> > extSpaceField;

    /** compute v_max, the maximum of the magnetization velocity of all the magnetic regions */
    void updateVmax(Eigen::VectorXd &sol);

    /** list of the Dirichlet indices where (vp,vq) are zero */
    std::vector<int> lvd;
    };// end class linAlgebra
#endif
