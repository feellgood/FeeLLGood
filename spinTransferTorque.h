#ifndef spinTransferTorque_h
#define spinTransferTorque_h

#include "config.h"
#include "pt3D.h"

/** \struct STT
 container for Spin Transfert Torque constants, Thiaville model, Dirichlet boundary conditions (potential fixed value on two or more surfaces)
 */

struct STT
    {
    double beta;/**< \f$ \beta \f$ is polarization rate of the current */    
    double N0;/**< density of states at Fermi level, units : J^-1 nm^-3  */
    double sigma;/**< Conductivity Ohm^-1 nm^-1 */
    double lJ;/**< length */
    double lsf;/**< spin flip length */
        
    std::vector<std::pair<std::string,double> > boundaryCond; /**< boundary conditions, first is the surface region name, second the associated value  */
    };
    
#endif
