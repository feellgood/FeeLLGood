#ifndef spinTransferTorque_h
#define spinTransferTorque_h

#include "config.h"
#include "pt3D.h"

/** \struct STT
 container for Spin Transfert Torque constants, Thiaville model, Dirichlet boundary conditions (potential fixed value on two or more surfaces)
 */

struct STT
    {
    int reg;/**< volume region number */
    double beta;/**< \f$ \beta \f$ is polarization rate of the current */    
    double N0;/**< density of states at Fermi level, units : J^-1 nm^-3  */
    double sigma;/**< Conductivity Ohm^-1 nm^-1 */
    double lJ;/**< length */
    double lsf;/**< spin flip length */
    std::function<double (Pt::pt3D)> func;/**< function to take into account spacial variation of current density (input is gauss point) */
        
    std::vector<std::pair<std::string,double> > boundaryCond; /**< boundary conditions, first is the surface region name, second the associated value  */
    
    /** electrostatic potential */ 
    std::vector<double> V; // size of the vector of nodes
    
    
    /** setter : node potential */
    inline void setNodesPotential(read_vector const& Xr)
        { 
        for (unsigned int i=0; i < V.size(); i++) V[i] = Xr[i];
        }
    
    /** getter for electrostatic potential */
    inline double get_elec_pot(const int i) const {return V[i];}
    
    /** setter for electrostatic potential */
    inline void set_elec_pot(const int i,const double val) {V[i] = val;}
    
    };
    
#endif
