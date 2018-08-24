#ifndef facette_h
#define facette_h

/** \file facette.h
  \brief contains namespace Facette
  header containing Fac class, and some constants and a less_than operator to redo orientation of triangular faces
 */

#include "feellgoodSettings.h"

#include "node.h"

/** \namespace Facette
 to grab altogether some constants for struct Fac
 */
namespace Facette
{
const int N = 3; /**< number of sommits */
const int NPI = 4; /**< number of weights  */

const double u[NPI]   = {   1/3.,   1/5.,   3/5.,   1/5.};/**< some constants to build hat functions */
const double v[NPI]   = {   1/3.,   1/5.,   1/5.,   3/5.};/**< some constants  to build hat functions */
const double pds[NPI] = {-27/96., 25/96., 25/96., 25/96.};/**< some constant weights  to build hat functions */




/** \class Fac
Face is a class containing the index references to nodes, it has a triangular shape and should not be degenerated 
*/
class Fac{
	public:
		inline Fac() {reg = 0;} /**< default constructor */
		int reg;/**< .msh region number */
		double surf; /**< surface of the face */
		double Ms; /**< magnetization at saturation of the face */    
		Pt::pt3D n;/**< normal vector n=(x,y,z) */	
		int ind[N];/**< indices table */
		double weight[NPI];/**< weights table */
		double a[N][NPI];          /**< hat functions table */
    
		/** computes the integral contribution of the triangular face */
		void integrales(Settings &mySets,std::vector <Node> &myNode, std::vector <double> &BE);
		
		/**
		convenient getter for N, usefull for templates projection and assemblage
		*/		
		inline int getN(void) {return N;}	
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
