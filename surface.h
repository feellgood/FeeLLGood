#ifndef surface_h
#define surface_h

/** \file surface.h
  \brief contains namespace Surface
  header containing Surf class
 */
 
 #include "facette.h"
 
 
 /** \namespace Facette
 to grab altogether some constants and calculation functions for class Surf
 */
namespace Surface
{

const int DIM = 2; /**< surface mathematical dimension */

class Surf{
	public:
		/** constructor used by readMesh */
        inline Surf(const std::vector<Nodes::Node> & _p_node /**< [in] vector of nodes */,
                   const int _reg /**< [in] region number */,
                   const int _idx /**< [in] region index in surface region vector */) : idxPrm(_idx),reg(_reg),refNode(_p_node) {}
	
		/** getter for region */
        	inline int getRegion(void) const {return reg;}
	
		/** getter for name */
        	inline std::string getName(void) const {return name;}
        	
        	/** push back method to fill the idx vector of facettes */
		inline void push_back( const int i ) { idx.push_back(i); }
	
	private:
		std::string name;/**< physical gmsh name */
		std::vector<int> idx;/**< indices of the facettes of the surface */
		
		const int idxPrm;/**< surface region index in surface region vector */
		
		const int reg;/**< .msh region number */
        
        	/** direct access to the Nodes */
        	const std::vector<Nodes::Node> & refNode;
		
};//end class Surf

}//end namespace
 
#endif /* surface_h */

