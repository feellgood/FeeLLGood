#ifndef facette_h
#define facette_h

/** \file facette.h
  \brief contains namespace Facette
  header containing Fac class, and some constants and a less_than operator to redo orientation of triangular faces
 */

#include "config.h"
#include "node.h"

/** \namespace Facette
 to grab altogether some constants and calculation functions for class Fac
 */
namespace Facette
{
const int N = 3; /**< number of sommits */
const int NPI = 4; /**< number of weights  */

const double u[NPI]   = {   1/3.,   1/5.,   3/5.,   1/5.};/**< some constants to build hat functions */
const double v[NPI]   = {   1/3.,   1/5.,   1/5.,   3/5.};/**< some constants  to build hat functions */
const double pds[NPI] = {-27/96., 25/96., 25/96., 25/96.};/**< some constant weights  to build hat functions */

/** hat function constants */
const double a[N][NPI] ={
    {1.-u[0]-v[0],1.-u[1]-v[1],1.-u[2]-v[2],1.-u[3]-v[3]},{u[0],u[1],u[2],u[3]},{v[0],v[1],v[2],v[3]}};

// for (int j=0; j<NPI; j++) { a[0][j] = 1.-u[j]-v[j]; a[1][j] = u[j]; a[2][j] = v[j];}


/** \class prm
region number and material constants
*/
struct prm
	{
	int reg;/**< region number */	
	double Js;/**< surface exchange */
	double Ks;/**< uniaxial surface anisotropy constant */	
	double uk[DIM]; /**< anisotropy axis */	
	
	/**
	 print the struct parameters
	 */
	inline void infos()
		{
		std::cout<< "surface region number = " << reg <<std::endl;
		std::cout<< "Js = " << Js <<std::endl;
		
		if(Ks!=0)
			{std::cout<< "Ks = " << Ks <<std::endl;
			std::cout<< "a = [ " << uk[0] << "," << uk[1] <<"," << uk[2] << "]" << std::endl;
			}
		else std::cout << "no surface anisotropy" << std::endl;
		};	
		
	};


/** \class Fac
Face is a class containing the index references to nodes, it has a triangular shape and should not be degenerated 
*/
class Fac{
	public:
		inline Fac() {reg = 0;} /**< default constructor */
		
		/** constructor from a region number and three indices */		
		inline Fac(int r,int i0,int i1,int i2) {reg = r; ind[0]=i0;ind[1]=i1;ind[2]=i2;}
		
		/** constructor from a region number, idxPrm and three indices */		
		inline Fac(int r,int idx,int i0,int i1,int i2) {reg = r; idxPrm=idx; ind[0]=i0;ind[1]=i1;ind[2]=i2;}
		
		int reg;/**< .msh region number */
		int idxPrm;/**< index of the material parameters of the facette */		
		double surf; /**< surface of the face */
		double Ms; /**< magnetization at saturation of the face */    
		Pt::pt3D n;/**< normal vector to the face */	
		int ind[N];/**< indices table of the nodes */
		double weight[NPI];/**< weights table */
		
        /** initialize weight  */
        inline void init(void)
            {for (int j=0; j<NPI; j++) {weight[j] = 2.*surf*pds[j]; }}// detJ = 2*surf;
    
		/** computes the integral contribution of the triangular face */
		void integrales(std::vector<Facette::prm> const& params,std::vector <Node> const& myNode, std::vector <double> &BE);
		
		/**
		convenient getter for N, usefull for templates projection and assemblage
		*/		
		inline int getN(void) {return N;}	

		/**
		computes normal vector n and surface surf		
		*/		
		void calc_surf(std::vector<Node> const& myNode);	
	};

/**
operator less_than for the orientation of the facette, lexicographic order
*/
struct less_than
{
/**
operator() for the comparison of two faces with lexicographical order
*/
bool operator()(Fac f1, Fac f2) const
  {
  if (f1.ind[0]<f2.ind[0]) return true;
  else
     if ((f1.ind[0]==f2.ind[0]) && (f1.ind[1]<f2.ind[1])) return true;
     else
        if ((f1.ind[0]==f2.ind[0]) && (f1.ind[1]==f2.ind[1]) && (f1.ind[2]<f2.ind[2])) return true;

  return false;
  }
};

}

#endif /* facette_h */
