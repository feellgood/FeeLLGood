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

/** \class LinAlgebra
convenient class to grab altogether some part of the calculations involved using eigen BiCGSTAB solver at each
timestep. The solver is handled by solver method, and is using Eigen::SparseMatrix, Row major matrix. This matrix
is prepared in 'batch mode', using a vector of triplets (also called COO write sparse matrix).
*/
class LinAlgebra
    {
public:
    /** constructor */
    inline LinAlgebra(Settings &s /**< [in] */, Mesh::mesh &my_msh /**< [in] */)
        : NOD(my_msh.getNbNodes()), MAXITER(s.MAXITER), TOL(s.TOL), ILU_tol(s.ILU_tol),
          ILU_fill_factor(s.ILU_fill_factor), verbose(s.verbose),
          prmTetra(s.paramTetra), prmFacette(s.paramFacette), refMsh(&my_msh)
        {
        Eigen::setNbThreads(s.solverNbTh);
        base_projection();
        if (!s.recenter)
            { idx_dir = Nodes::IDX_UNDEF; }
        else
            { idx_dir = s.recentering_direction; }

        if(s.getFieldType() == R4toR3)
            { setExtSpaceField(s); }
        }

    /** computes inner data structures of tetraedrons and triangular facettes (K matrices and L vectors) */
    void prepareElements(Eigen::Vector3d const &Hext /**< [in] applied field */, timing const &t_prm /**< [in] */);

    /** computes inner data structures of tetraedrons and triangular facettes (K matrices and L vectors) */
    void prepareElements(double const A_Hext /**< [in] amplitude applied field */, timing const &t_prm /**< [in] */);

    /** build init guess for bicgstab solver */
    void buildInitGuess(Eigen::Ref<Eigen::VectorXd> G) const;

    /** solver, uses eigen stabilized biconjugate gradient solver (bicgstab) with ILU preconditionner, sparse matrix and vector are filled    with multiThreading. Sparse matrix is row major.
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

    /** maximum number of iteration for bicgstab */
    const int MAXITER;

    /** solver tolerance */
    const double TOL;

    /** ILU preconditionner tolerance */
    double ILU_tol;

    /** ILU preconditionner filling factor */
    double ILU_fill_factor;

    /** verbosity */
    const int verbose;

    /** material parameters of the tetrahedrons */
    const std::vector<Tetra::prm> &prmTetra;

    /** material parameters of the facettes */
    const std::vector<Facette::prm> &prmFacette;

    /** direct access to the mesh */
    Mesh::mesh *refMsh;

    /** speed of the domain wall */
    double DW_vz;

    /** maximum speed of the magnetization in the whole physical object */
    double v_max;

    /** external applied space field, values on gauss points, size is number of tetraedrons */
    std::vector< Eigen::Matrix<double,Nodes::DIM,Tetra::NPI> > extSpaceField;
    };// end class linAlgebra
#endif
