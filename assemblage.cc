#include "fem.h"

void assemblage(Fem &fem, Tet &elt,
           gmm::dense_matrix <double> &Ke, vector <double> &Le,
           write_matrix &K, write_vector &L)
{
const int NOD = fem.NOD;
const int N   = Tet::N;

for (int i=0; i<N; i++){
    int i_= elt.ind[i];             
    for (int j=0; j<N; j++){
        int j_= elt.ind[j];
        K(NOD+i_,j_)+= Ke(i,j);      K(NOD+i_, NOD+j_)+= Ke(  i,N+j);
        K(    i_,j_)+= Ke(N+i,j);    K(    i_, NOD+j_)+= Ke(N+i,N+j);
	    }
    L[NOD+i_]+= Le[  i];
    L[    i_]+= Le[N+i];
    }
//cout<<"temps assemblage: "<<diff_t<<endl;
}

void assemblage(Fem &fem, Fac &elt,
           gmm::dense_matrix <double> &Ke, vector <double> &Le,
           write_matrix &K, write_vector &L)
{
const int NOD = fem.NOD;
const int N   = Fac::N;

for (int i=0; i<N; i++){
    int i_= elt.ind[i];             
    for (int j=0; j<N; j++){
        int j_= elt.ind[j];
        K(NOD+i_,j_)+= Ke(i,j);      K(NOD+i_, NOD+j_)+= Ke(  i,N+j);
        K(    i_,j_)+= Ke(N+i,j);    K(    i_, NOD+j_)+= Ke(N+i,N+j);
	    }
    L[NOD+i_]+= Le[  i];
    L[    i_]+= Le[N+i];
    }
//cout<<"temps assemblage: "<<diff_t<<endl;
}
