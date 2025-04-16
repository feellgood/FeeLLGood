#ifndef spinAccumulationSolver_h
#define spinAccumulationSolver_h

#include <vector>
#include "config.h"
#include "node.h"
#include "tetra.h"
#include "electrostatSolver.h"

/** \class spinAcc
 container for Spin Accumulation constants, diffusive accumulation spin model, Boundary conditions
 (potential fixed value on one surface, current density on another surface)
 */

class spinAcc
    {
    /** constructor */
    spinAcc(Mesh::mesh &_msh /**< [in] ref to the mesh */,
    std::vector<Tetra::prm> _pTetra /**< [in] ref to vector of param tetra (volume region parameters) */):
        msh(_msh), paramTetra(_pTetra)
        {
        }

    /** mesh */
    Mesh::mesh msh;

    private:
    /** this vector contains the material parameters for all regions for all the tetrahedrons */
    std::vector<Tetra::prm> paramTetra;

    /** \f$ \beta \f$ is polarization rate of the current */
    double getBeta(Tetra::Tet &tet) const;

    /** density of states at Fermi level, units : J^-1 nm^-3  */
    double getN0(Tetra::Tet &tet) const;

    /** length */
    double getLsd(Tetra::Tet &tet) const;

    /** spin flip length */
    double getLsf(Tetra::Tet &tet) const;

    /** affect extraField function and extraCoeffs_BE function for all the tetrahedrons */
    void prepareExtras(std::vector<Tetra::Tet> &v_tet, electrostatSolver &elec);
    };

#endif
