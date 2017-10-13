#include "fem.h"

void init_distrib(Fem &fem) 
{
const int NOD = fem.NOD;
//const double L=220.;
//const double scale = fem.scale;

for (int i=0; i<NOD; i++) {
   Node &node = fem.node[i];
   //double x = node.x/scale;
   //double y = node.y/scale;
   //double z = node.z/scale;
// double r = sqrt(x*x+y*y+z*z+1e-8);
   
   node.u[0] = 1./sqrt(2.);
                    // fabs(cos(2.*M_PI/L*x))*sin(0.01);
   node.u[1] = 0.;
                    // sin(2.*M_PI/L*x);
   node.u[2] = 1./sqrt(2.);
                    // fabs(cos(2.*M_PI/L*x))*cos(0.01);
   node.phi  = 0.;
   }     
}
