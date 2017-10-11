#include "fem.h"

void projection(Fem &fem, Tet &elt,
           gmm::dense_matrix <double> &A,  vector <double> &B,
           gmm::dense_matrix <double> &Ap, vector <double> &Bp)
{
const int N = Tet::N;
gmm::dense_matrix <double> P(2*N,3*N), PA(2*N,3*N);
for (int i=0; i<N; i++){
    Node &node = fem.node[elt.ind[i]];
    P(i,i)  = node.ep[0];  P(i,N+i)  = node.ep[1];  P(i,2*N+i)  = node.ep[2];
    P(N+i,i)= node.eq[0];  P(N+i,N+i)= node.eq[1];  P(N+i,2*N+i)= node.eq[2];
    }

mult(P,A,PA);
mult(PA, gmm::transposed(P), Ap);

mult(P,B,Bp);
}

void projection(Fem &fem, Fac &elt,
           gmm::dense_matrix <double> &A,  vector <double> &B,
           gmm::dense_matrix <double> &Ap, vector <double> &Bp)
{
const int N = Fac::N;
gmm::dense_matrix <double> P(2*N,3*N), PA(2*N,3*N);
for (int i=0; i<N; i++){
    Node &node = fem.node[elt.ind[i]];
    P(i,i)  = node.ep[0];  P(i,N+i)  = node.ep[1];  P(i,2*N+i)  = node.ep[2];
    P(N+i,i)= node.eq[0];  P(N+i,N+i)= node.eq[1];  P(N+i,2*N+i)= node.eq[2];
    }

mult(P,A,PA);
mult(PA, gmm::transposed(P), Ap);

mult(P,B,Bp);
}
