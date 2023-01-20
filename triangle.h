#ifndef triangle_h
#define triangle_h

/** \file triangle.h
  \brief contains Triangle class
  Triangle class is a basic object containing three indices refering to nodes
 */

#include "node.h"

namespace Mesh
{
const int TRIANGLE_NB_SOMMITS = 3; /**< number of sommits */

/** \class Triangle
geometrical triangle : three indices refering to nodes
*/
class Triangle
    {
    public:
		/** constructor used by readMesh */
        inline Triangle(const std::vector<Nodes::Node> & _p_node /**< [in] vector of nodes */,
                   const int i0 /**< [in] node index */,
                   const int i1 /**< [in] node index */,
                   const int i2 /**< [in] node index */) : refNode(_p_node)
        {
        ind[0]=i0;
        ind[1]=i1;
        ind[2]=i2;
        
        int _NOD = _p_node.size();
        
        if (_NOD>0)
            { // to force index to start from 0 (C++) instead of Matlab/msh convention
            if((i0>0)&&(i0<=_NOD)) 
            	{ ind[0]--; } 
            else 
            	{ std::cout<< "warning index i0 out of bounds in Triangle constructor" <<std::endl; }
            	
            if((i1>0)&&(i1<=_NOD))
            	{ ind[1]--; }
            else 
            	{ std::cout<< "warning index i1 out of bounds in Triangle constructor" <<std::endl; }
            
            if((i2>0)&&(i2<=_NOD))
            	{ ind[2]--; }
            else
            	{ std::cout<< "warning index i2 out of bounds in Triangle constructor" <<std::endl; }
            }
        }

		int ind[TRIANGLE_NB_SOMMITS];/**< indices table of the nodes */
	
	private:
		/** direct access to the Nodes */
        	const std::vector<Nodes::Node> & refNode;
	};

} //end namespace
#endif
