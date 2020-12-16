#ifndef mesh_h
#define mesh_h

#include <set>

#include "node.h"
#include "facette.h"
#include "tetra.h"

#include "feellgoodSettings.h"

class mesh
{
public:
    /** default constructor */
    inline mesh() {}
    
    /** constructor */
    inline mesh(Settings & mySets /**< [in] */)
        {
        readMesh(mySets);
        femutil(mySets);// reordering of index nodes for facette orientation, also some modifications on fac::Ms
        geometry();// initialization of l,c,diam,vol,surf
            //init_distrib(mySets);
        }
    
    inline int getNbNodes(void) { return node.size(); }
    inline int getNbFacs(void) { return fac.size(); }
    inline int getNbTets(void) { return tet.size(); }
    
    inline const Nodes::Node getNode(const int i) { return node[i]; }
    
    inline Nodes::Node & setNode(const int i) {return node[i];}
    
    /** basic informations on the mesh */
    void infos(void) const
        {
        std::cout << "\t diam bounding box ="<< diam << std::endl;
        std::cout << "\t nodes\t\t\t" << node.size() << std::endl;
        std::cout << "\t faces\t\t\t" << fac.size() << std::endl;
        std::cout << "\t tetraedrons\t\t" << tet.size() << std::endl;
        std::cout << "\t Total surface\t\t"  << surf << std::endl;
        std::cout << "\t Total volume\t\t\t" << vol << std::endl;
        }
    
    
    /** check phi values for debug */
    inline void checkNodes(void)
        {
        int k=0;
        for(unsigned int i=0;i<node.size();i++) { if (std::isnan(node[i].phi)) {std::cout << "#" << i<<std::endl; node[i].phi = 0;k++;}  }
        std::cout << "total nb of NaN in node[].phi ="<< k << std::endl;
        }

    inline void evolution(void) { std::for_each(node.begin(), node.end(), [](Nodes::Node &n){ n.evolution();} ); }
    
    /** find center and length along coordinates and diameter = max(l(x|y|z)), computes vol,surf */
    void geometry(void)
        {
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
        }
        
    Pt::pt3D c;/**< center position */	
	Pt::pt3D l;/**< lengths along x,y,z axis */	
	
	double diam;/**< max of l coordinates, to define a bounding box */
	double surf;/**< total surface */
	double vol;/**< total volume of the mesh */
    
    std::vector <Facette::Fac>  fac; /**< face container */
	std::vector <Tetra::Tet>  tet; /**< tetrahedron container */
    
    /**
    reset the nodes struct to restart another step time simulation
    Undo the action of one or many "vsolve" runs in case of failure.
    Demagnetizing field and energies don't need to be reset, because they won't be updated if failure is detected.
    I don't know how to cleanly reset "fem.DW_vz". BC
    */
    inline void reset(void) { std::for_each(node.begin(),node.end(),[](Nodes::Node &n) {n.reset();}); }
    
    
    
    /**
read a solution from a file (tsv formated) and initialize fem struct to restart computation from that distribution, return time
*/
double readSol(bool VERBOSE/**< [in] */,
             double scaling /**< [in] scaling factor for physical coordinates */,
             const std::string fileName /**< [in] */ );

/** computes an analytical initial magnetization distribution as a starting point for the simulation */
    inline void init_distrib(Settings & mySets /**< [in] */)
        { std::for_each( node.begin(),node.end(), [this,&mySets](Nodes::Node & n) 
            {
            Pt::pt3D pNorm = Pt::pt3D( (n.p.x() - c.x())/l.x() , (n.p.y() - c.y())/l.y() , (n.p.z() - c.z())/l.z() );
            n.u0 = mySets.getValue(pNorm);
            n.u = n.u0; n.phi  = 0.; n.phiv = 0.;
            } ); 
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
    
    void calc_charges(std::function<const Pt::pt3D (Nodes::Node)> getter,std::vector<double> & srcDen,std::vector<double> & corr,Settings const& settings)
    {
    int nsrc = 0;
    std::for_each(tet.begin(),tet.end(),[&srcDen,getter,&nsrc,&settings](Tetra::Tet const& tet)              
        { tet.charges(getter,srcDen,nsrc, nu0 * settings.paramTetra[tet.idxPrm].J ); });//end for_each on tet

    std::for_each(fac.begin(),fac.end(),[&srcDen,&corr,getter,&nsrc](Facette::Fac const& fac)
        { fac.charges(getter,srcDen,corr,nsrc); });// end for_each on fac
    }
    
private:
    std::vector <Nodes::Node> node; /**< node container */
	
	/** reading mesh file function */
    void readMesh(Settings const& mySets);
	
	/** read old mesh format 2.2 */
    void readOldMesh(Settings const& mySets,std::ifstream &msh);

    /** read mesh format 4.1 */
    void readNewMesh(Settings const& mySets,std::ifstream &msh);
	
	/** return the minimum of all nodes coordinate along coord axis */
    inline double minNodes(const Pt::index coord)
        {
        const auto minCoord = std::min_element(node.begin(),node.end(),[coord](Nodes::Node &n1,Nodes::Node &n2) {return (n1.p(coord)<n2.p(coord)); } );
        return minCoord->p(coord); 
        }

    /** return the maximum of all nodes coordinate along coord axis */
    inline double maxNodes(const Pt::index coord)
        {
        const auto maxCoord = std::max_element(node.begin(),node.end(),[coord](Nodes::Node &n1,Nodes::Node &n2) {return (n1.p(coord)<n2.p(coord)); } );
        return maxCoord->p(coord);
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
    void femutil(Settings const& settings)
        {
        std::set<Facette::Fac> sf;//implicit use of operator< overloaded in class Fac

        std::for_each(tet.begin(),tet.end(),[&sf](Tetra::Tet const& te)
            {
            int ia=te.ind[0];int ib=te.ind[1];int ic=te.ind[2];int id=te.ind[3];
            const int reg = te.getRegion();
            sf.insert( Facette::Fac(reg,te.idxPrm,ia,ic,ib) );
            sf.insert( Facette::Fac(reg,te.idxPrm,ib,ic,id) );
            sf.insert( Facette::Fac(reg,te.idxPrm,ia,id,ic) );
            sf.insert( Facette::Fac(reg,te.idxPrm,ia,ib,id) ); 
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
                        Facette::Fac fc(node.size());
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
