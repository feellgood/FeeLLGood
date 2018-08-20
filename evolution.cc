#include "fem.h"

#include <algorithm>

void Fem::evolution(void)
{
/*
//const int NOD = fem.NOD;
for (int i=0; i<NOD; i++) {
    //Node &node = fem.node[i];
    for (int d=0; d<3; d++) {
        node[i].u0[d] = node[i].u[d];
        node[i].v0[d] = node[i].v[d];
    }
    node[i].phi0  = node[i].phi;
    node[i].phiv0 = node[i].phiv;
}
*/

std::for_each(node.begin(), node.end(), [](Node &n){ for(int d=0;d<D;d++) {n.u0[d]=n.u[d];n.v0[d]=n.v[d];}; n.phi0=n.phi; n.phiv0=n.phiv; } );


for (int e=0; e<4; e++) //fem.E0[e] = fem.E[e];
	E0[e] = E[e];
Etot0 = Etot;//fem.Etot0 = fem.Etot;
}
