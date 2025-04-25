#ifndef spinAccumulationSolver_h
#define spinAccumulationSolver_h

#include <vector>
#include "config.h"

/** \struct STT
 container for Spin Transfert Torque constants, Thiaville model, Dirichlet boundary conditions
 (potential fixed value on two or more surfaces)
 */

struct STT
    {
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

    /** if true a text file V.sol is generated. It contains the solution of the
                     electrostatic problem on the nodes of the mesh */
    bool V_file;

    /** boundary conditions, stored as a vector of pairs.
    First element of the pair is the surface region name given in the mesh by its physical name;
    Second element is the electrostatic potential associated value  */
    std::vector<std::pair<std::string, double>> boundaryCond;
    };

#endif
