#ifndef spinAccumulationSolver_h
#define spinAccumulationSolver_h

#include <vector>
#include "config.h"
#include "node.h"
#include "tetra.h"
#include "electrostatSolver.h"

/** \class spinAcc
 container for Spin Accumulation constants, diffusive accumulation spin model, Dirichlet boundary conditions
 (potential fixed value on two or more surfaces)
 */

class spinAcc
    {
    /** constructor */
    spinAcc(Mesh::mesh &_msh, double beta, double N0, double _sigma, double lJ, double lsf):
        msh(_msh) //, elec(_msh,1e-6,false,5000,false,"")
        {
        ksi = Nodes::sq(lJ / lsf);
        D0 = 2.0 * _sigma / (Nodes::sq(CHARGE_ELECTRON) * N0);
        pf = Nodes::sq(lJ) / (D0 * (1. + ksi * ksi)) * BOHRS_MUB * beta / CHARGE_ELECTRON;
        }

    /** mesh */
    Mesh::mesh msh;

    /** \f$ \beta \f$ is polarization rate of the current */
    double beta;

    /** density of states at Fermi level, units : J^-1 nm^-3  */
    double N0;

    /** Conductivity Ohm^-1 nm^-1 */
    double sigma;

    /** length */
    double lJ;

    /** spin flip length */
    double lsf;

    /** solve electrostatic side problem */
    //electrostatSolver elec;

    inline double calc_prefactor(double Js) { return D0 / Nodes::sq(lJ) / (gamma0*Js/mu0); }

    /** ksi is in Thiaville notations beta_DW */
    double ksi;

    /** density of states */
    double D0;

    /** a prefactor for BE coefficient coefficients*/
    double pf;

    /** affect extraField function and extraCoeffs_BE function for all the tetrahedrons */
    void prepareExtras(std::vector<Tetra::Tet> &v_tet, electrostatSolver &elec);
    };

#endif
