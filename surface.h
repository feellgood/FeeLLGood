#ifndef surface_h
#define surface_h

/** \file surface.h
  \brief contains namespace Surface
  header containing Surf class
 */
 
 #include "triangle.h"
 
 
namespace Mesh
{

const int SURF_DIM = 2; /**< surface mathematical dimension */


/** \class Surf
Surf is a class containing a vector of triangular facettes, it has to be oriented
*/
class Surf{
	public:
		/** constructor used by readMesh */
        	inline Surf(const std::vector<Nodes::Node> & _p_node /**< [in] vector of nodes */,
                   const int _reg /**< [in] region number from mesh */,
                   const std::string _name /**< [in] physical name from mesh */
                   ) : name(_name),reg(_reg),refNode(_p_node) {}
	
		/** getter for region */
        	inline int getRegion(void) const {return reg;}
	
		/** getter for name */
        	inline std::string getName(void) const {return name;}
        	
        	/** push back method to fill the vector of triangles */
		inline void push_back( Mesh::Triangle t ) { elem.push_back(t); }
		
		/** getter for the number of triangular facettes */
		inline int getNbElem(void) { return elem.size(); }
	
	private:
		const std::string name;/**< physical gmsh name */
		
		const int reg;/**< .msh region number */
        
        	/** direct access to the Nodes */
        	const std::vector<Nodes::Node> & refNode;
		
		std::vector<Mesh::Triangle> elem;/**< triangles of the surface */
		
};//end class Surf

}//end namespace
 
#endif /* surface_h */

