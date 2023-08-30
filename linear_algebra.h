#ifndef linear_algebra_h
#define linear_algebra_h

/** \file linear_algebra.h
\brief secondary header, it grabs altogether the linear algebra by the solver to apply fem method
<br> It encapsulates the calls to GMM, the assemblage and projection of the matrix for all elements
<br> projection and matrix assembly is multithreaded for tetrahedron, monothread for facette
*/
#include <random>
#include <execution>

#include "config.h"

#include "facette.h"
#include "feellgoodSettings.h"
#include "mesh.h"
#include "node.h"
#include "pt3D.h"
#include "tetra.h"

/** \class LinAlgebra
convenient class to grab altogether some part of the calculations involved using gmm solver at each
timestep
*/
class LinAlgebra
    {
public:
    /** constructor */
    inline LinAlgebra(Settings &s /**< [in] */, Mesh::mesh &my_msh /**< [in] */)
        : NOD(my_msh.getNbNodes()), refMsh(&my_msh), settings(s)
        {
        base_projection();
        }

    /** computes inner data structures of tetraedrons and triangular facettes (K matrices and L
     * vectors) */
    void prepareElements(Pt::pt3D const &Hext /**< [in] applied field */, timing const &t_prm /**< [in] */);

    /**  solver, uses bicgstab and gmres, sparse matrix and vector are filled with multiThreading */
    int solver(timing const &t_prm /**< [in] */);

    /** setter for DW_dz */
    inline void set_DW_vz(double vz /**< [in] */) { DW_vz = vz; }

    /** getter for v_max */
    inline double get_v_max(void) { return v_max; }

private:
    const int NOD; /**< total number of nodes, also an offset for filling sparseMatrix, initialized
                      by constructor */
    Mesh::mesh *refMsh; /**< direct access to the mesh */

    double DW_vz;             /**< speed of the domain wall */
    const Settings &settings; /**< settings */
    double v_max;             /**< maximum speed */

    /** computes local vector basis {ep,eq} in the tangeant plane for projection on the elements */
    void base_projection();

    };  // fin class linAlgebra

#endif
