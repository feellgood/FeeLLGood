/** \file linear_algebra.h
\brief secondary header, it grabs altogether the linear algebra by the solver to apply fem method  <br>
It encapsulates the calls to GMM , the assemblage and projection of the matrix for all elements <br>
two templates projection and assemblage template class parameter is either Facette::Fac or Tetra::Tet 
*/

#include "boost/numeric/mtl/mtl.hpp"
#include "boost/numeric/itl/itl.hpp"

#include <boost/numeric/mtl/matrix/inserter.hpp>
#include <boost/numeric/mtl/operation/set_to_zero.hpp>
#include <boost/numeric/mtl/interface/vpt.hpp>

#include "fem.h"

#ifndef linear_algebra_h
#define linear_algebra_h

/** convenient typedef for mtl4 */
typedef typename mtl::Collection< mtl::compressed2D<double> >::value_type v_type;

/** convenient typedef for mtl4 inserter */
typedef mtl::mat::inserter< mtl::compressed2D<double>,mtl::update_plus<v_type> > sparseInserter;


/** \class LinAlgebra
convenient class to grab altogether some part of the calculations involved using gmm solver at each timestep
*/
class LinAlgebra
{
public:
	/** constructor */	
	inline LinAlgebra(Settings & s,
                      std::vector<Node> & myNode,
                      std::vector <Tetra::Tet> & myTet,
                      std::vector <Facette::Fac> & myFace) 
    {settings = &s; refNode = &myNode; refTet = &myTet; refFac = &myFace;}
	
	/** pointer to diagonal preconditionner  */
	itl::pc::diagonal < mtl::compressed2D<double> > *prc;

    /** solver, uses bicgstab and gmres */
	int  vsolve(double dt,long nt);

    /** setter for DW_dz */
    inline void set_DW_vz(double vz){DW_vz = vz;}    

    /** setter for Hext */
    inline void set_Hext(double Hx,double Hy,double Hz){Hext[0]=Hx;Hext[1]=Hy;Hext[2]=Hz;}

    /** getter for v_max */
    inline double get_v_max() {return v_max;}
    
    /** getter node */
    inline Node getNode(int i) {return (*refNode)[i];}
    
    /** getter node physical position */
    inline Pt::pt3D getNodePhysPos(int i) {return (*refNode)[i].p;} 
    
private:
    std::vector<Node>  *refNode;/**< direct access to the Nodes */
	std::vector <Facette::Fac> *refFac; /**< direct access to the faces */
	std::vector <Tetra::Tet> *refTet; /**< direct access to the tetrahedrons */
	
	double Hext[DIM];/**< applied field */
    double DW_vz;/**< speed of the domain wall */
	Settings *settings;/**< pointer to the settings */
    double v_max;/**< maximum speed */
	
/** computes the local vector basis {ep,eq} in the tangeant plane for projection on the elements */
inline void base_projection(void)
	{ std::for_each(refNode->begin(),refNode->end(),[](Node &n) { n.buildBase_epeq();}); }
	

    void feedMat(const int NOD,double dt, mtl::compressed2D<double> &K_T,mtl::dense_vector<double> &L_T,std::vector<Tetra::Tet>::iterator it_b,std::vector<Tetra::Tet>::iterator it_e);
    
}; // fin class linAlgebra

#endif
