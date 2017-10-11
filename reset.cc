#include "fem.h"

/*
 * Undo the action of one or many "vsolve" runs in case of failure.
 * Demagnetizing field and energies don't need to be reset, because they won't
 * be updated if failure is detected.
 * I don't know how to cleanly reset "fem.DW_vz".
 */

void reset(Fem &fem)
{
const int NOD = fem.NOD;

for (int i=0; i<NOD; i++) {
    Node &node = fem.node[i];
    for (int d=0; d<3; d++) {
        node.u[d] = node.u0[d];
        node.v[d] = node.v0[d];
    }
    node.phi  = node.phi0;
    node.phiv = node.phiv0;
}

}
