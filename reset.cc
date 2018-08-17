#include "fem.h"

/*
 * Undo the action of one or many "vsolve" runs in case of failure.
 * Demagnetizing field and energies don't need to be reset, because they won't
 * be updated if failure is detected.
 * I don't know how to cleanly reset "fem.DW_vz".
 */

void Fem::reset(void)
{
//const int NOD = fem.NOD;

for (int i=0; i<NOD; i++) {
    //Node &node = fem.node[i];
    for (int d=0; d<D; d++) {
        node[i].u[d] = node[i].u0[d];
        node[i].v[d] = node[i].v0[d];
    }
    node[i].phi  = node[i].phi0;
    node[i].phiv = node[i].phiv0;
}

}
