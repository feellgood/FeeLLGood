/** \file fem.h
\brief principal header, contains the struct fem <br>
This file is called by almost all source files in feellgood project <br>
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

#include <time.h>
#include <sys/times.h>
#include <unistd.h>

#include "config.h"
#include "feellgoodSettings.h"

#include "ANN.h" // ANN declarations

#include "pt3D.h"
#include "tetra.h"
#include "facette.h" // tiny.h est includé là

#include "node.h"

const bool U = true;/**< used as a template parameter */
const bool V = false;/**< used as a template parameter */


/** \struct Fem
container to grab altogether all parameters of a simulation, including mesh geometry, containers for the mesh
*/
struct Fem{
	int NOD;/**< number of nodes in vector container */
	//int SEQ;/**< number of sequences, usefull to define a vector applied field eventually varying in time */
	Pt::pt3D c;/**< center position */	
	Pt::pt3D l;/**< lengths along x,y,z axis */	
	
	double diam;/**< max of l coordinates, to define a bounding box */
	double surf;/**< total surface */
	double vol;/**< total volume of the mesh */
	
	double fmm_normalizer;/**< normalizing geometrical constant for the computation of the demag field */
	
	double t;/**< physical current time of the simulation */
	
	double vmax;/**< maximum speed of magnetization */

	double E0[4];/**< table to store initial energy values <br> 
index convention : 0-exchange 1-anisotropy 2-demagnetizing 3-applied */

	double E[4]; /**< table to store energy values at time t <br>
index convention : 0-exchange 1-anisotropy 2-demagnetizing 3-applied */

	double DW_vz0;/**< initial speed of the domain wall (to check) */
	double DW_vz;/**< speed of the domain wall along z */
	double DW_dir;/**< direction of the domain wall (to check) */
	double DW_z; /**< domain wall displacement along Oz */
	
	double Etot0;/**< initial total energy (to check) */
	double Etot;/**< total energy */
	double evol;/**< increment dE for dt */
	double phy;/**< no idea */
	std::vector <Nodes::Node> node; /**< node container */
	std::vector <Facette::Fac>  fac; /**< face container */
	std::vector <Tetra::Tet>  tet; /**< tetrahedron container */
    
    triple Hext;/**< external applied field direction (should be normalized) */
    
    ANNkd_tree* kdtree;/**< ANN kdtree to find efficiently the closest set of nodes to a physical point in the mesh  */
    ANNpointArray pts;/**< container for the building of the kdtree (handled by ANN library) */
    

/**
print some informations of fem container
*/
inline void infos(void)
{
std::cout << "This is feeLLGood SHA1= " + std::string(SHAnumber) << std::endl;
std::cout << "diam bounding box ="<< diam << std::endl;
std::cout << "\t nodes\t\t\t" << node.size() << std::endl;
std::cout << "\t faces\t\t\t" << fac.size() << std::endl;
std::cout << "\t tetraedrons\t\t" << tet.size() << std::endl;
std::cout << "\t Total surface\t\t"  << surf << std::endl;
std::cout << "\t Total volume\t\t\t" << vol << std::endl;
}

/**
reset the nodes struct to restart another step time simulation
Undo the action of one or many "vsolve" runs in case of failure.
Demagnetizing field and energies don't need to be reset, because they won't be updated if failure is detected.
I don't know how to cleanly reset "fem.DW_vz". BC
*/
inline void reset(void) { std::for_each(node.begin(),node.end(),[](Nodes::Node &n) {n.reset();}); }

/**
time evolution : one step in time
*/
inline void evolution(void)
{
std::for_each(node.begin(), node.end(), [](Nodes::Node &n){ n.evolution();} );
for (int e=0; e<4; e++) { E0[e] = E[e]; }
Etot0 = Etot;
}

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
inline void init_distrib(Settings &mySets)
	{ std::for_each( node.begin(),node.end(), [this,&mySets](Nodes::Node &n) 
        {
        Pt::pt3D pNorm = Pt::pt3D( (n.p.x() - c.x())/l.x() , (n.p.y() - c.y())/l.y() , (n.p.z() - c.z())/l.z() );
        n.u = mySets.getValue(pNorm);
        n.phi  = 0.;} 
    ); }

/**
read a solution from a file (tsv formated) and initialize fem struct to restart computation from that distribution
*/
void readSol(double scaling /**< [in] scaling factor for physical coordinates */,
	std::string fileName /**< [in] */ );

/** reading mesh file function */
void readMesh(Settings &mySets);

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
void femutil(Settings &settings /**< [in] */);

/** find direction of motion of DW */
void direction(enum Pt::index idx_dir);

/** computes energies stored in E table */
void energy(Settings &settings);

/** recentering algorithm for the study of the motion of a micromagnetic object (domain wall). 
 
if \f$ D_i>0 \f$				        if  \f$ D_i<0 \f$

<----------------|------->		------->|<----------------	\f$ m_i = < u_i > < 0 \f$

or					or

---------------->|<-------		<-------|---------------->	\f$ m_i = <u_i> > 0 \f$
*/
bool recentrage(double thres/**< [in] threshold parameter */,enum Pt::index idx_dir /**< [in] */);

/** saving function for a solution */
void saver(Settings &settings, std::ofstream &fout, int nt);

/** text file (vtk) writing function for a solution */
void savecfg_vtk(std::string fileName);

/** text file (tsv) writing function for a solution */
void savesol(std::string fileName,double s);

/** save the field values */
void saveH(std::string fileName,double scale);


/** 
template to compute average component of either u or v on the whole set of tetetrahedron
*/
template <int UorV>
double moy(Pt::index d)
{
double sum = 0.;
std::for_each(tet.begin(),tet.end(),[this,&sum,&d](Tetra::Tet &te)
    {
    double val_nod[Tetra::N], val[Tetra::NPI];
    for (int ie=0; ie<Tetra::N; ie++) 
        {
        int i = te.ind[ie];
        Nodes::Node &n = node[i];
        if(UorV) { val_nod[ie] = n.u(d);} else { val_nod[ie] = n.v(d);} 
        }
    tiny::transposed_mult<double, Tetra::N, Tetra::NPI> (val_nod, Tetra::a, val);
    //sum += tiny::sp<double, Tetra::NPI> (val, te.weight);
    sum += te.weightedScalarProd(val);
    }
);//fin for_each

return sum/vol;
}


}; // end struct fem definition
