#include <algorithm>

#include "linear_algebra.h"


/*----------------------------------------------------------------------*/
/* Construction de la base de projection du plan tangent  	    	*/
/*----------------------------------------------------------------------*/

void LinAlgebra::base_projection(void) 
{
std::for_each(fem.node.begin(),fem.node.end(),
[](Node &n) {
    Pt::pt3D u = n.u0;
    //double u0 = n.u0[0];
    //double u1 = n.u0[1];
    //double u2 = n.u0[2]; 
    //double norme2 = Pt::pScal(u,u);//u0*u0 + u1*u1 + u2*u2;
// cout<< "nod " << i << ": u0=" << u0<< ", u1=" << u1<< ", u2=" << u2<< endl;
    if ( Pt::norme2(u) > 0){            // vrai si u est defini       
        Pt::pt3D r = Pt::rand();
        //double r1 = rand() / (RAND_MAX+1.);
        //double r2 = rand() / (RAND_MAX+1.);  
//cout << "r0=" << r0 << ", r1=" << r1 << ", r2=" << r2 << endl;
        //double t0 = r1*u2-r2*u1;
        //double t1 = r2*u0-r0*u2;
        //double t2 = r0*u1-r1*u0;
	Pt::pt3D t = r*u;        
	t.normalize();	
	//double norme = sqrt(t0*t0 + t1*t1 + t2*t2);
        //t0/=norme; t1/=norme; t2/=norme;

// premier vecteur tangent       
	n.ep = t;        
	//n.ep[0] = t0; n.ep[1] = t1; n.ep[2] = t2;
//cout<< "nod " << i << ": t0=" << t0<< ", t1=" << t1<< ", t2=" << t2< endl;
// second vecteur tangent              
	n.eq = u*t;        
	//n.eq[0] = u1*t2 - u2*t1; n.eq[1] = u2*t0 - u0*t2; n.eq[2] = u0*t1 - u1*t0;

//cout << "nod " << i << ": eq0=" << node.eq.x << ", eq1=" << node.eq.y <<
//        ", eq2=" << node.eq.z << endl;
        }
    else{ n.ep = Pt::pt3D(0.,0.,0.); n.eq = Pt::pt3D(0.,0.,0.);} }
);// fin du for_each
    
}
