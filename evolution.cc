#include "fem.h"

void evolution(Fem &fem)
{
const int NOD = fem.NOD;

for (int i=0; i<NOD; i++) {
    Node &node = fem.node[i];
    for (int d=0; d<3; d++) {
        node.u0[d] = node.u[d];
        node.v0[d] = node.v[d];
    }
    node.phi0  = node.phi;
    node.phiv0 = node.phiv;
}

for (int e=0; e<4; e++)
    fem.E0[e] = fem.E[e];
fem.Etot0 = fem.Etot;
}
