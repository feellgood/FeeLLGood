#ifndef mesh_h
#define mesh_h

/** \file mesh.h
\brief class mesh, read_mesh, save
*/

#include <set>
#include <execution>
#include <map>

#include "node.h"
#include "facette.h"
#include "tetra.h"
#include "surface.h"

#include "feellgoodSettings.h"


namespace Mesh {

/** \class mesh
class for storing the mesh, including mesh geometry values, containers for the nodes, triangular faces and tetrahedrons.
 nodes data are private. They are accessible only through getter and setter.
*/
class mesh
{
public:
    /** constructor : read mesh file, reorder indices and computes some values related to the mesh :
     center and length along coordinates and diameter = max(l(x|y|z)), volume and surface */
    inline mesh(Settings const& mySets /**< [in] */)
        {
        readMesh(mySets);
        indexReorder(mySets);// reordering of index nodes for facette orientation, also some modifications on fac::Ms
        
        double xmin = minNodes(Pt::IDX_X);
        double xmax = maxNodes(Pt::IDX_X);

        double ymin = minNodes(Pt::IDX_Y);
        double ymax = maxNodes(Pt::IDX_Y);

        double zmin = minNodes(Pt::IDX_Z);
        double zmax = maxNodes(Pt::IDX_Z);

        l = Pt::pt3D(xmax - xmin,ymax - ymin,zmax - zmin);
        diam = l.maxLength();
        c = Pt::pt3D(0.5*(xmax + xmin),0.5*(ymax + ymin),0.5*(zmax + zmin));
        
        vol = std::transform_reduce(std::execution::par,tet.begin(),tet.end(),0.0, std::plus{},
                                    [](Tetra::Tet const& te) {return te.calc_vol();} );

        surf = std::transform_reduce(std::execution::par, fac.begin(), fac.end(),0.0, std::plus{},
                                    [](Facette::Fac const& fa) {return fa.calc_surf();} );
        }
    
    /** return number of nodes  */
    inline int getNbNodes(void) const { return node.size(); }
    
    /** return number of triangular fac */
    inline int getNbFacs(void) const { return fac.size(); }
    
    /** return number of tetrahedrons */
    inline int getNbTets(void) const { return tet.size(); }
    
    /** getter : return node */
    inline const Nodes::Node getNode(const int i) const { return node[i]; }
    
    /** setter for u0 */
    inline void set_node_u0(const int i,Pt::pt3D const& val) { node[i].u0 = val; }
    
    /** fix to zero node[i].v */
    inline void set_node_zero_v(const int i) { node[i].v = Pt::pt3D(0.,0.,0.); }

    /** basic informations on the mesh */
    void infos(void) const
        {
        std::cout << "mesh:\n";
        std::cout << "  bounding box diam:  " << diam << '\n';
        std::cout << "  nodes:              " << getNbNodes() << '\n';
        std::cout << "  faces:              " << fac.size() << '\n';
        std::cout << "  tetraedrons:        " << tet.size() << '\n';
        std::cout << "  total surface:      " << surf << '\n';
        std::cout << "  total volume:       " << vol << '\n';
        }

/** call setBasis for all nodes */
void setBasis(const double r)
	{
	std::for_each(std::execution::par,node.begin(),node.end(), [&r](Nodes::Node &nod){ nod.setBasis(r); } );
	}

/** make_evol on all nodes, and returns v_max */
double updateNodes(std::vector<double> const& X,const double dt)
	{
	double v2max = 0.0;
	const unsigned int NOD = node.size();

	for(unsigned int i=0; i < NOD ; i++)
	    {
	    double vp = X[i]*gamma0;
	    double vq = X[NOD+i]*gamma0;
	    double v2 = vp*vp + vq*vq;
	    if (v2>v2max) { v2max = v2; }
	    node[i].make_evol(vp,vq,dt);    
	    }

	return sqrt(v2max);
	}


    /** call evolution for all the nodes */
   inline void evolution(void) 
   {
   std::for_each(std::execution::par,node.begin(),node.end(),[](Nodes::Node &nod) { nod.evolution(); } );
   }
    
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
    
    /** surface container */
	std::vector <Mesh::Surf>  s;
    

    /** read a solution from a file (tsv formated) and initialize fem struct to restart computation from that distribution, return time
    */
    double readSol(bool VERBOSE/**< [in] */, const std::string fileName /**< [in] input .sol text file */ );

    /** computes an analytical initial magnetization distribution as a starting point for the simulation */
    inline void init_distrib(Settings const& mySets /**< [in] */)
        { 
        for(unsigned int i=0;i< node.size();i++)
            {
            node[i].u0 = mySets.getValue( Nodes::get_p(node[i]) );
            node[i].u = node[i].u0; node[i].phi  = 0.; node[i].phiv = 0.;
            }  
        }
    
    /** 
    average component of either u or v through getter on the whole set of tetetrahedron
    */
    double avg(std::function<double (Nodes::Node,Pt::index)> getter /**< [in] */,Pt::index d /**< [in] */) const
    {
    double sum = std::transform_reduce(std::execution::par,tet.begin(),tet.end(),0.0, std::plus{}, [getter,&d](Tetra::Tet const& te)
    	{
    	double val[Tetra::NPI]; 
        te.interpolation(getter,d,val); 
        return te.weightedScalarProd(val);
    	});
    
    return sum/vol;
    }

    /** text file (vtk) writing function for a solution, using text VTK format (unstructured grid version 2.0) : deprecated, it is recommanded to use convert2vtk python script instead */
    void savecfg_vtk(timing const& t_prm /**< [in] */,const std::string fileName /**< [in] */) const;

    /** text file (tsv) writing function for a solution */
    void savesol(const int precision /**< [in] numeric precision in .sol output text file */,
    const std::string fileName /**< [in] */,std::string const& metadata /**< [in] */) const;

    /** text file (tsv) writing function for a solution of a side problem, used by electrostatSolver */
    bool savesol(const int precision /**< [in] */, const std::string fileName /**< [in] */,
     std::string const& metadata /**< [in] */, std::vector<double> const& val /**< [in] */) const;

    /** save the demagnetizing field values, including idx and npi indices, for debug use */
    void saveH(const std::string fileName /**< [in] */,const double t/**< [in] */,const double scale /**< [in] */) const;
    
    /** computes all charges for the demag field to feed a tree in the fast multipole algo (scalfmm) */
    void calc_charges(std::function<const Pt::pt3D (Nodes::Node)> getter,std::vector<double> & srcDen,std::vector<double> & corr,Settings const& settings)
    {
    int nsrc = 0;
    std::for_each(tet.begin(),tet.end(),[&srcDen,getter,&nsrc,&settings](Tetra::Tet const& tet)              
        { tet.charges(getter,srcDen,nsrc, nu0 * settings.paramTetra[tet.idxPrm].J ); });

    std::for_each(fac.begin(),fac.end(),[&srcDen,&corr,getter,&nsrc](Facette::Fac const& fac)
        { fac.charges(getter,srcDen,corr,nsrc); });
    }
    
    /** setter for node[i]; what_to_set will fix what is the part of the node struct to set (usefull for fmm_demag.h) */
    inline void set(const int i, std::function<void (Nodes::Node &,const double)> what_to_set,const double val)
    	{ what_to_set(node[i],val); }
    
private:

    /** node container: not initialized by constructor, but later while reading the mesh by member function init_node */
    std::vector< Nodes::Node > node;
    
    /** map of the surface region physical names from mesh file */
    std::map< int, std::string > surfRegNames;
    
    /** map of the volume region physical names from mesh file */
    std::map< int, std::string > volRegNames;
    
    /** memory allocation for the nodes */
    inline void init_node(const int Nb) { node.resize(Nb); }
    
	/** reading mesh format 2.2 text file function */
    void readMesh(Settings const& mySets);

    /** loop on nodes to apply predicate 'whatTodo'  */
    double doOnNodes(const double init_val,const Pt::index coord, std::function<bool(double,double)> whatToDo) const
        {
        double result(init_val);
        std::for_each(node.begin(),node.end(), [&result,coord,whatToDo](Nodes::Node const& n)
                    {
                    double val = n.p(coord);
                    if( whatToDo(val,result)) result = val;
                    }  );
        return result;
        }
    
	/** return the minimum of all nodes coordinate along coord axis */
    inline double minNodes(const Pt::index coord) const
        {
        return doOnNodes(__DBL_MAX__ , coord, [](double a, double b){return a<b;} );
        }

    /** return the maximum of all nodes coordinate along coord axis */
    inline double maxNodes(const Pt::index coord) const
        {
        return doOnNodes(__DBL_MIN__ , coord, [](double a, double b){return a>b;} );
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
                        Facette::Fac fc(node,0,0,0,0,0);
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
        if(settings.verbose)
            { std::cout << "  reindexed\n"; }
        }
};

} //end namespace Mesh

#endif
