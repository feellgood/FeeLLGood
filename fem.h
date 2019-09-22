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

#include "config.h"
#include "feellgoodSettings.h"

#include "ANN.h" // ANN declarations

#include "pt3D.h"
#include "tetra.h"
#include "facette.h"

#include "node.h"

const bool U = true;/**< used as a template parameter */
const bool V = false;/**< used as a template parameter */

/** \class Fem
class container to grab altogether all parameters of a simulation, including mesh geometry, containers for the mesh
*/
class Fem
    {
    public:
        /** constructor */
        inline Fem(Settings & mySets) //mySets cannot be passed const because of getValue method in init_distrib, to set x,y,z values
            {
            readMesh(mySets);
            femutil(mySets);
            chapeaux(mySets.EPSILON);
            
            if (mySets.restore)
                { readSol(mySets.verbose,mySets.getScale(), mySets.restoreFileName); }
            else
                {
                std::cout<< "initial magnetization(x,y,z,t=0) :\nMx =" << mySets.sMx << "\nMy =" << mySets.sMy << "\nMz =" << mySets.sMz << std::endl; 
                init_distrib(mySets);
                }    
            direction(mySets.verbose,Pt::IDX_Z);/* determination de la direction de propagation de la paroi */
            t=0.;
            
            Hext[0] = nu0*mySets.Bext[0];
            Hext[1] = nu0*mySets.Bext[1];
            Hext[2] = nu0*mySets.Bext[2]; 
            
            vmax  = 0.0;
            DW_vz = DW_vz0 = 0.0;
            DW_z  = 0.0;
            fmm_normalizer = 1./(2.*diam);
            E_exch0 = E_exch = 0.0;
            E_aniso0 = E_aniso = 0.0;
            E_demag0 = E_demag = 0.0;
            E_zeeman0 = E_zeeman = 0.0;
            Etot0 = Etot = 0.0;
            }
        
	Pt::pt3D c;/**< center position */	
	Pt::pt3D l;/**< lengths along x,y,z axis */	
	
	double diam;/**< max of l coordinates, to define a bounding box */
	double surf;/**< total surface */
	double vol;/**< total volume of the mesh */
	
	double fmm_normalizer;/**< normalizing geometrical constant for the computation of the demag field */
	
	double t;/**< physical current time of the simulation */
	
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
	double evol;/**< increment dE for dt */
	
	std::vector <Nodes::Node> node; /**< node container */
	std::vector <Facette::Fac>  fac; /**< face container */
	std::vector <Tetra::Tet>  tet; /**< tetrahedron container */
    
    triple Hext;/**<  field. Unit : SI (A/m) */
    
    /**
    print some informations of fem container
    */
    void infos(void) const;
    
    /** computes all the energies */
    void energy(Settings const& settings /**< [in] */);
    
    /**
    time evolution : one step in time
    */
    inline void evolution(void)
        {
        std::for_each(node.begin(), node.end(), [](Nodes::Node &n){ n.evolution();} );
        E_exch0 = E_exch;
        E_aniso0 = E_aniso;
        E_demag0 = E_demag;
        E_zeeman0 = E_zeeman;
        Etot0 = Etot;
        }
        
    /**
    reset the nodes struct to restart another step time simulation
    Undo the action of one or many "vsolve" runs in case of failure.
    Demagnetizing field and energies don't need to be reset, because they won't be updated if failure is detected.
    I don't know how to cleanly reset "fem.DW_vz". BC
    */
    inline void reset(void) { std::for_each(node.begin(),node.end(),[](Nodes::Node &n) {n.reset();}); }    
    
    /** saving function for a solution */
    void saver(Settings const& settings /**< [in] */, std::ofstream &fout /**< [out] */, const int nt /**< [in] */) const;

    /** text file (vtk) writing function for a solution */
    void savecfg_vtk(Settings const& settings /**< [in] */,const std::string fileName /**< [in] */) const;

    /** text file (tsv) writing function for a solution */
    void savesol(const std::string fileName /**< [in] */,const double s /**< [in] */) const;

    /** save the field values */
    void saveH(const std::string fileName /**< [in] */,const double scale /**< [in] */) const;
    
    /** 
    average component of either u or v through getter on the whole set of tetetrahedron
    */
    double avg(std::function<double (Nodes::Node,Pt::index)> getter /**< [in] */,Pt::index d /**< [in] */) const
    {// syntaxe pénible avec opérateur binaire dans la lambda pour avoir un += sur la fonction voulue, with C++17 we should use reduce instead of accumulate here
    double sum = std::accumulate(tet.begin(),tet.end(),0.0, [&getter,&d](double &s,Tetra::Tet const& te)
                            {
                            double val[Tetra::NPI]; 
                            te.interpolation(getter,d,val); 
                            return (s + te.weightedScalarProd(val));    
                            } );

    return sum/vol;
    }
    
    /** recentering algorithm for the study of the motion of a micromagnetic object (domain wall). 
 
    if \f$ D_i>0 \f$				        if  \f$ D_i<0 \f$

    <----------------|------->		------->|<----------------	\f$ m_i = < u_i > < 0 \f$

    or					or

    ---------------->|<-------		<-------|---------------->	\f$ m_i = <u_i> > 0 \f$
    */
    bool recentrage(double thres/**< [in] threshold parameter */,enum Pt::index idx_dir /**< [in] */);

    
    
    private:
    ANNkd_tree* kdtree;/**< ANN kdtree to find efficiently the closest set of nodes to a physical point in the mesh  */
    ANNpointArray pts;/**< container for the building of the kdtree (handled by ANN library) */
    
/**
computes the hat functions for all containers
*/
inline void chapeaux(double epsilon /**< [in] if \f$ | detJ | < \epsilon \f$ jacobian is considered degenerated */) 
{
std::for_each(fac.begin(),fac.end(),[](Facette::Fac &f) {f.init();});

std::for_each(tet.begin(),tet.end(),[epsilon](Tetra::Tet &t) { t.init(epsilon); });
}
/**
computes an analytical initial magnetization distribution as a starting point for the simulation
*/
inline void init_distrib(Settings & mySets /**< [in] */)
	{ std::for_each( node.begin(),node.end(), [this,&mySets](Nodes::Node & n) 
        {
        Pt::pt3D pNorm = Pt::pt3D( (n.p.x() - c.x())/l.x() , (n.p.y() - c.y())/l.y() , (n.p.z() - c.z())/l.z() );
        n.u0 = mySets.getValue(pNorm);// u or u0?
        n.u = n.u0;
        n.phi  = 0.;} 
    ); }

/**
read a solution from a file (tsv formated) and initialize fem struct to restart computation from that distribution
*/
void readSol(bool VERBOSE/**< [in] */,
             double scaling /**< [in] scaling factor for physical coordinates */,
             std::string fileName /**< [in] */ );

/** read old mesh format 2.2 */
void readOldMesh(Settings const& mySets,std::ifstream &msh);

/** read old mesh format 4.1 */
void readNewMesh(Settings const& mySets,std::ifstream &msh);

/** reading mesh file function */
void readMesh(Settings const& mySets);

/**  utilitary function to initialize pts,kdtree,l,c,diam, computes the surfaces and volumes and reorientation of the tetrahedrons if needed in fem struct; definition of Ms on facette elements <br>
Indices and orientation convention : 

                        v
                      .
                    ,/
                   /
                2(ic)                                 2
              ,/|`\                                 ,/|`\
            ,/  |  `\                             ,/  |  `\
          ,/    '.   `\                         ,6    '.   `5
        ,/       |     `\                     ,/       8     `\
      ,/         |       `\                 ,/         |       `\
     0(ia)-------'.--------1(ib) --> u     0--------4--'.--------1
      `\.         |      ,/                 `\.         |      ,/
         `\.      |    ,/                      `\.      |    ,9
            `\.   '. ,/                           `7.   '. ,/
               `\. |/                                `\. |/
                  `3(id)                                `3
                     `\.
                        ` w

*/
        void femutil(Settings const& settings /**< [in] */);

        /** find direction of motion of DW */
        void direction(bool VERBOSE /**< [in] VERBOSE mode */, enum Pt::index idx_dir /**< [in] */);


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
