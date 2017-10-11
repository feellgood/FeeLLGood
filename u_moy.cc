#include "fem.h"
#include "tiny.h"

/*-----------------------------------------*/
/*   valeur moyenne d'une composante       */
/*-----------------------------------------*/

double u_moy(Fem &fem, int d)
{
const int TET = fem.TET;
const int N   = Tet::N;
const int NPI = Tet::NPI;

double sum = 0.;
for (int t=0; t<TET; t++){
    Tet &tet = fem.tet[t];
    double val_nod[N], val[NPI];
    for (int ie=0; ie<N; ie++) {
        int i = tet.ind[ie];
        Node &node = fem.node[i];
        val_nod[ie] = node.u[d];
        }
   tiny::transposed_mult<double, N, NPI> (val_nod, tet.a, val);
   sum += tiny::sp<double, NPI> (val, tet.weight);
   }

return sum/fem.vol;
}
