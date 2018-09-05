#include <algorithm>

#include "linear_algebra.h"


/*----------------------------------------------------------------------*/
/* Construction de la base de projection du plan tangent  	    	*/
/*----------------------------------------------------------------------*/

void LinAlgebra::base_projection(void) 
{
std::for_each(fem.node.begin(),fem.node.end(),[](Node &n) { n.buildBase_epeq();});// fin du for_each
    
}
