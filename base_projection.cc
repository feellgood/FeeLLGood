#include "fem.h"

/*----------------------------------------------------------------------*/
/* Construction de la base de projection du plan tangent  	    	*/
/*----------------------------------------------------------------------*/

void base_projection(Fem &fem) 
{
const int NOD = fem.NOD;
    
for (int i=0; i<NOD; i++){
    Node &node = fem.node[i];
    double u0,u1,u2, t0,t1,t2;
    u0 = node.u0[0];
    u1 = node.u0[1];
    u2 = node.u0[2]; 
    double norme2 = u0*u0 + u1*u1 + u2*u2;
// cout<< "nod " << i << ": u0=" << u0<< ", u1=" << u1<< ", u2=" << u2<< endl;
    if (norme2 > 0){            // vrai si u est defini       
        double r0 = rand() / (RAND_MAX+1.);
        double r1 = rand() / (RAND_MAX+1.);
        double r2 = rand() / (RAND_MAX+1.);  
//cout << "r0=" << r0 << ", r1=" << r1 << ", r2=" << r2 << endl;
        t0 = r1*u2-r2*u1;
        t1 = r2*u0-r0*u2;
        t2 = r0*u1-r1*u0;
        double norme = sqrt(t0*t0 + t1*t1 + t2*t2);
        t0/=norme; t1/=norme; t2/=norme;

// premier vecteur tangent       
        node.ep[0] = t0;
        node.ep[1] = t1;
        node.ep[2] = t2;
//cout<< "nod " << i << ": t0=" << t0<< ", t1=" << t1<< ", t2=" << t2< endl;
// second vecteur tangent              
        node.eq[0] = u1*t2 - u2*t1;
        node.eq[1] = u2*t0 - u0*t2;
        node.eq[2] = u0*t1 - u1*t0;

//cout << "nod " << i << ": eq0=" << node.eq.x << ", eq1=" << node.eq.y <<
//        ", eq2=" << node.eq.z << endl;
        }
    else{
        node.ep[0] = node.ep[1] = node.ep[2] = 0;
        node.eq[0] = node.eq[1] = node.eq[2] = 0;}}
}
