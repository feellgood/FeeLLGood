#ifndef mesh_h
#define mesh_h

/** \file mesh.h
\brief class mesh, read_mesh, save
*/

#include <set>

#include "node.h"
#include "facette.h"
#include "tetra.h"

#include "feellgoodSettings.h"

/** \class mesh
class for storing the mesh, including mesh geometry values, containers for the nodes, triangular faces and tetrahedrons. nodes data are not public. They are accessible only through getter and setter.
*/
class mesh
{
public:
    /** constructor  : read mesh, reorder indices and computes values related to the mesh : center and length along coordinates and diameter = max(l(x|y|z)), vol,surf */
    inline mesh(Settings const& mySets /**< [in] */):nbNod(0)
        {
        readMesh(mySets);
        if(mySets.verbose) std::cout<< "mesh in memory." <<std::endl;
        indexReorder(mySets);// reordering of index nodes for facette orientation, also some modifications on fac::Ms
        if(mySets.verbose) std::cout<< "mesh reindexed." <<std::endl;
        
        double xmin = minNodes(Pt::IDX_X);
        double xmax = maxNodes(Pt::IDX_X);

        double ymin = minNodes(Pt::IDX_Y);
        double ymax = maxNodes(Pt::IDX_Y);

        double zmin = minNodes(Pt::IDX_Z);
        double zmax = maxNodes(Pt::IDX_Z);

        l = Pt::pt3D(xmax - xmin,ymax - ymin,zmax - zmin);
        diam = l.maxLength();
        c = Pt::pt3D(0.5*(xmax + xmin),0.5*(ymax + ymin),0.5*(zmax + zmin));
        vol = std::accumulate(tet.begin(),tet.end(),0.0,[](double x,Tetra::Tet const& te){return x + te.calc_vol();} );
        surf = std::accumulate(fac.begin(),fac.end(),0.0,[](double x,Facette::Fac const& fa){return x + fa.calc_surf();} );
        
        if(mySets.verbose) std::cout<< "mesh geometry computed." <<std::endl;
        }
    
    /** return number of nodes  */
    inline int getNbNodes(void) const { return nbNod; }
    
    /** return number of triangular fac */
    inline int getNbFacs(void) const { return fac.size(); }
    
    /** return number of tetrahedrons */
    inline int getNbTets(void) const { return tet.size(); }
    
    /** getter : return node */
    inline const Nodes::Node getNode(const int i) const { return node[i]; }
    
    /** setter : node */
    inline Nodes::Node & setNode(const int i) {return node[i];}
    
    /** setter : node potential */
    inline void setNodesPotential(read_vector const& Xr)
        { for (int i=0; i < nbNod; i++) node[i].V = Xr[i]; }
    
    /** basic informations on the mesh */
    void infos(void) const
        {
        std::cout << "\t diam bounding box ="<< diam << std::endl;
        std::cout << "\t nodes\t\t\t" << getNbNodes() << std::endl;
        std::cout << "\t faces\t\t\t" << fac.size() << std::endl;
        std::cout << "\t tetraedrons\t\t" << tet.size() << std::endl;
        std::cout << "\t Total surface\t\t"  << surf << std::endl;
        std::cout << "\t Total volume\t\t\t" << vol << std::endl;
        }

    /** call evolution for all the nodes */
   inline void evolution(void) { for(int i=0;i<nbNod;i++) node[i].evolution(); }
    
    /** isobarycenter */	
    Pt::pt3D c;

    /** lengths along x,y,z axis */
    Pt::pt3D l;
	
    /** max of l coordinates, to define a bounding box */
	double diam;
    
    /** total surface */
	double surf;
    
    /** total volume of the mesh */
	double vol;
    
     /** face container */
    std::vector <Facette::Fac>  fac;
    
     /** tetrahedron container */
	std::vector <Tetra::Tet>  tet;
    
    /**
    reset the nodes struct to restart another step time simulation
    Undo the action of one or many "vsolve" runs in case of failure.
    Demagnetizing field and energies don't need to be reset, because they won't be updated if failure is detected.
    I don't know how to cleanly reset "fem.DW_vz". BC
    */
    inline void reset(void) { for(int i=0;i<nbNod;i++) node[i].reset(); }
    

    /** read a solution from a file (tsv formated) and initialize fem struct to restart computation from that distribution, return time
    */
    double readSol(bool VERBOSE/**< [in] */,
             double scaling /**< [in] scaling factor for physical coordinates */,
             const std::string fileName /**< [in] */ );

    /** computes an analytical initial magnetization distribution as a starting point for the simulation */
    inline void init_distrib(Settings const& mySets /**< [in] */)
        { 
        for(int i=0;i<nbNod;i++)
            {
            node[i].u0 = mySets.getValue( Nodes::get_p(node[i]) );
            node[i].u = node[i].u0; node[i].phi  = 0.; node[i].phiv = 0.;
            }  
        }
    
    /** 
    average component of either u or v through getter on the whole set of tetetrahedron
    */
    double avg(std::function<double (Nodes::Node,Pt::index)> getter /**< [in] */,Pt::index d /**< [in] */) const
    {// syntaxe pénible avec opérateur binaire dans la lambda pour avoir un += sur la fonction voulue, with C++17 we should use reduce instead of accumulate here
    double sum = std::accumulate(tet.begin(),tet.end(),0.0, [getter,&d](double &s,Tetra::Tet const& te)
                            {
                            double val[Tetra::NPI]; 
                            te.interpolation(getter,d,val); 
                            return (s + te.weightedScalarProd(val));    
                            } );
    return sum/vol;
    }

    /** text file (vtk) writing function for a solution, using text VTK format (unstructured grid version 2.0) : deprecated, it is recommanded to use convert2vtk python script instead */
    void savecfg_vtk(Settings const& settings /**< [in] */,timing const& t_prm /**< [in] */,const std::string fileName /**< [in] */) const;

    /** text file (tsv) writing function for a solution */
    void savesol(const std::string fileName /**< [in] */,timing const& t_prm /**< [in] */,const double s /**< [in] */) const;

    /** save the demagnetizing field values, including idx and npi indices, for debug use */
    void saveH(const std::string fileName /**< [in] */,const double t/**< [in] */,const double scale /**< [in] */) const;
    
    /** computes all charges for the demag field to feed a tree in the fast multipole algo (scalfmm) */
    void calc_charges(std::function<const Pt::pt3D (Nodes::Node)> getter,std::vector<double> & srcDen,std::vector<double> & corr,Settings const& settings)
    {
    int nsrc = 0;
    std::for_each(tet.begin(),tet.end(),[&srcDen,getter,&nsrc,&settings](Tetra::Tet const& tet)              
        { tet.charges(getter,srcDen,nsrc, nu0 * settings.paramTetra[tet.idxPrm].J ); });//end for_each on tet

    std::for_each(fac.begin(),fac.end(),[&srcDen,&corr,getter,&nsrc](Facette::Fac const& fac)
        { fac.charges(getter,srcDen,corr,nsrc); });// end for_each on fac
    }
    
    /** getter for electrostatic potential */
    inline double get_elec_pot(const int i) const {return node[i].V;}
    
    /** setter for electrostatic potential */
    inline void set_elec_pot(const int i,const double val) {node[i].V = val;}
    
private:
    /** node container : the node list is not initialized by constructor, but later, while reading the mesh, by member function init_node */
    Nodes::NodeList node;
    
    
    
    /** memory allocation for the nodes */
    inline void init_node(const int Nb)
        {
        nbNod = Nb;
        node = Nodes::NodeList(new Nodes::Node[Nb]);
        }
    
    /** total number of nodes read from mesh file */
    int nbNod;
    
    
	/** reading mesh format 2.2 file function */
    void readMesh(Settings const& mySets);

    
	/** return the minimum of all nodes coordinate along coord axis */
    inline double minNodes(const Pt::index coord) const
        {
        double _min = __DBL_MAX__;
        for(int i=0;i<nbNod;i++)
            {
            double val = node[i].p(coord);
            if (val < _min) _min = val;
            }
        return _min;
        }

    /** return the maximum of all nodes coordinate along coord axis */
    inline double maxNodes(const Pt::index coord) const
        {
        double _max = __DBL_MIN__;
        for(int i=0;i<nbNod;i++)
            {
            double val = node[i].p(coord);
            if (val > _max) _max = val;
            }
        return _max;
        }
    
    /** redefine orientation of triangular faces in accordance with the tetrahedron
    * reorientation of the tetrahedrons if needed; definition of Ms on facette elements
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
    void indexReorder(Settings const& settings)
        {
        std::set<Facette::Fac> sf;//implicit use of operator< overloaded in class Fac

        std::for_each(tet.begin(),tet.end(),[&sf](Tetra::Tet const& te)
            {
            std::set<Facette::Fac> tet_set = te.ownedFac();
            sf.insert(tet_set.begin(),tet_set.end());
            });//end for_each

        std::for_each(fac.begin(),fac.end(),[this,&settings,&sf](Facette::Fac &fa)
            {
            fa.Ms = 0.;
            if ( !(settings.paramFacette[fa.idxPrm].suppress_charges) )
                {
                int i0 = fa.ind[0], i1 = fa.ind[1], i2 = fa.ind[2];
                std::set< Facette::Fac>::iterator it=sf.end();
                for (int perm=0; perm<2; perm++) 
                    {
                    for (int nrot=0; nrot<3; nrot++)
                        {
                        Facette::Fac fc(getNbNodes());
                        fc.ind[(0+nrot)%3]=i0; fc.ind[(1+nrot)%3]=i1; fc.ind[(2+nrot)%3]=i2;
                        it=sf.find(fc);
                        if (it!=sf.end()) break;
                        }
      
                    if (it!=sf.end()) 
                        { // found
                        const Pt::pt3D & p0 = node[ it->ind[0] ].p;
                        const Pt::pt3D & p1 = node[ it->ind[1] ].p;
                        const Pt::pt3D & p2 = node[ it->ind[2] ].p;
                
                        //fa.Ms will have the magnitude of first arg of copysign, with the sign of second arg
                        fa.Ms = std::copysign(nu0*settings.paramTetra[it->idxPrm].J , Pt::pTriple( p1-p0 , p2-p0 ,fa.calc_norm()) );// carefull here, calc_norm recomputes the normal to the face before the idx swap
                        }
                    std::swap(i1,i2);// it seems from ref archive we do not want to swap inner fac indices but local i1 and i2
                    }//end perm
                }
            }); //end for_each
        }
};

#endif
