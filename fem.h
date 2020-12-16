#ifndef fem_h
#define fem_h

/** \file fem.h
\brief principal header, contains the struct fem
This file is called by many source files in feellgood.
It does also contains the definition of many constants for the solver, and for scalfmm 
*/
#include <iomanip>
#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <map>
#include <utility>
#include <cmath>
#include <cstdlib>
#include <cassert>
#include <list>
#include <functional>

#include <time.h>
#include <sys/times.h>
#include <unistd.h>

#include "mesh.h"

#include "feellgoodSettings.h"
#include "ANN.h"

/** \class Fem
class container to grab altogether all parameters of a simulation, including mesh geometry, containers for the mesh
*/
class Fem
    {
    public:
        /** constructor: call mesh constructor, initialize pts,kdtree and many inner variables */
        inline Fem(Settings & mySets,timing &t_prm):msh(mySets) //mySets cannot be passed const because of getValue method in init_distrib, to set x,y,z values
            {
            vmax  = 0.0;
            DW_vz = DW_vz0 = 0.0;
            DW_z  = 0.0;
            
            E_exch0 = E_exch = 0.0;
            E_aniso0 = E_aniso = 0.0;
            E_demag0 = E_demag = 0.0;
            E_zeeman0 = E_zeeman = 0.0;
            Etot0 = Etot = 0.0;
            
            pts= annAllocPts(msh.getNbNodes(), Pt::DIM);
            kdtree = new ANNkd_tree(pts, msh.getNbNodes(), Pt::DIM);
            if (!kdtree) SYSTEM_ERROR;
            
            for(int i=0;i<msh.getNbNodes();i++)
                { 
                Nodes::Node const& n = msh.getNode(i);
                this->pts[i][0] = n.p.x();
                this->pts[i][1] = n.p.y();
                this->pts[i][2] = n.p.z(); 
                }
            
            if (mySets.restoreFileName == "")
                {
                if(mySets.verbose)
                    { std::cout<< "initial magnetization M(x,y,z,t=0) = { " << mySets.sMx << "\t" << mySets.sMy << "\t" << mySets.sMz << " }\n"; } 
                msh.init_distrib(mySets);
                }
            else
                {
                double _t = msh.readSol(mySets.verbose,mySets.getScale(), mySets.restoreFileName);
                t_prm.set_t(_t);    
                }
                    
            direction(Pt::IDX_Z);/* determination de la direction de propagation de la paroi */
            }
    
    /** destructor */
    ~Fem ()
        { annDeallocPts(pts); delete kdtree; }
    
    
    double vmax;/**< maximum speed of magnetization */

	double E_exch0; /**< previous iteration exchange energy  */
    double E_aniso0; /**< previous iteration anisotropy energy  */
	double E_demag0; /**< previous iteration demagnetizing energy  */
    double E_zeeman0; /**< previous iteration zeeman energy  */

	double E_exch; /**< exchange energy */
    double E_aniso; /**< anisotropy energy  */
	double E_demag; /**< demagnetizing energy  */
    double E_zeeman; /**< zeeman energy  */
    
	double DW_vz0;/**< initial speed of the domain wall */
	double DW_vz;/**< speed of the domain wall along z */
	double DW_dir;/**< direction of the domain wall */
	double DW_z; /**< domain wall displacement along Oz */
	
	double Etot0;/**< initial total energy */
	double Etot;/**< total energy */
	
	mesh msh;
    
    /**
    print some informations of fem container
    */
    void infos(void) const;
    
    /** computes all the energies */
    void energy(double const t /**< [in] time in second, used to compute zeeman contribution if applied field is time dependant */,
                Settings & settings /**< [in] */);
    
    /**
    time evolution : one step in time
    */
    inline void evolution(void)
        {
        msh.evolution();
        E_exch0 = E_exch;
        E_aniso0 = E_aniso;
        E_demag0 = E_demag;
        E_zeeman0 = E_zeeman;
        Etot0 = Etot;
        }
        
    
    /** saving function for a solution */
    void saver(Settings & settings /**< [in] */,timing const& t_prm /**< [in] */, std::ofstream &fout /**< [out] */, const int nt /**< [in] */) const;
    
    /** recentering algorithm for the study of the motion of a domain wall. 
 
    if \f$ D_i>0 \f$				        if  \f$ D_i<0 \f$

    <----------------|------->		------->|<----------------	\f$ m_i = < u_i > < 0 \f$

    or					or

    ---------------->|<-------		<-------|---------------->	\f$ m_i = <u_i> > 0 \f$
    */
    
    
    bool recenter(double thres/**< [in] threshold parameter */,char recentering_direction /**< [in] X|Y|Z */);

    private:
    ANNkd_tree* kdtree;/**< ANN kdtree to find efficiently the closest set of nodes to a physical point in the mesh  */
    ANNpointArray pts;/**< container for the building of the kdtree (handled by ANN library) */

/** find direction of motion of DW */
void direction( enum Pt::index idx_dir /**< [in] */);


/** zeroing of all energies */
inline void zeroEnergy(void) {
        E_exch = 0.0;
        E_aniso = 0.0;
        E_demag = 0.0;
        E_zeeman = 0.0;
        Etot = 0.0;
        }

/** computes the sum of all energies */
inline void calc_Etot(void) { Etot = E_exch + E_aniso + E_demag + E_zeeman; }

}; // end class

#endif

