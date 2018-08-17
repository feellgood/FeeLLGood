#include "linear_algebra.h"

using namespace std;

void chapeaux(Fem &fem)
{
const int FAC = fem.FAC;
const int TET = fem.TET;

/********************* FACES *******************/
for (int f=0; f<FAC; f++){
   Fac &fac = fem.fac[f];
   
const int NPI = Fac::NPI;

// NPI 4
   double u[NPI]   = {   1/3.,   1/5.,   3/5.,   1/5.};
   double v[NPI]   = {   1/3.,   1/5.,   1/5.,   3/5.};
   double pds[NPI] = {-27/96., 25/96., 25/96., 25/96.};
  
   double detJ = 2*fac.surf;

   for (int j=0; j<NPI; j++){
       fac.a[0][j] = 1.-u[j]-v[j];
       fac.a[1][j] = u[j];
       fac.a[2][j] = v[j];
       fac.weight[j]    = detJ * pds[j];
       }
   }

/****************** TETRAS *****************/
for (int t=0; t<TET; t++){
    Tet &tet = fem.tet[t];
    const int N = Tet::N;
    const int NPI  = Tet::NPI;

// NPI 5
    double A,B,C,D,E;
    A=1./4.;  B=1./6.;  C=1./2.;  D=-2./15.;  E=3./40.;
    double u[NPI]   = {A,B,B,B,C};
    double v[NPI]   = {A,B,B,C,B};
    double w[NPI]   = {A,B,C,B,B};
    double pds[NPI] = {D,E,E,E,E}; 
    
    gmm::dense_matrix <double> nod(3,N); 
    for (int i=0; i<N; i++){
        int i_= tet.ind[i];
//	cout << "t:" << t <<",  nod:" << i_ << endl;
        nod(0,i) = fem.node[i_].x;
        nod(1,i) = fem.node[i_].y;
        nod(2,i) = fem.node[i_].z;
        }
// cout << nod <<endl;

    for (int j=0; j<NPI; j++){
        vector <double> a(N);
        gmm::dense_matrix <double> da(N,3), dadu(N,3), J(3,3);
        a[0] = 1.-u[j]-v[j]-w[j];
        a[1] = u[j];
        a[2] = v[j];
        a[3] = w[j];

        dadu(0,0)= -1;   dadu(0,1)= -1;   dadu(0,2)= -1;
        dadu(1,0)= +1;   dadu(1,1)=  0;   dadu(1,2)=  0;
        dadu(2,0)=  0;   dadu(2,1)= +1;   dadu(2,2)=  0;
        dadu(3,0)=  0;   dadu(3,1)=  0;   dadu(3,2)= +1;
	
        mult(nod, dadu, J);
        double detJ = lu_det(J);
//	cout << "tet: " << t << "   jac:" << detJ <<endl;
	
        if (fabs(detJ) < EPSILON){
#ifdef LIBRARY
            ostringstream what;
            what << "Singular jacobian in tetrahedron " << t;
            throw runtime_error(what.str());
#else
            cerr << "jacobienne singuliere ds le tetraedre " << t << endl;
            exit(1);
#endif
            }

        lu_inverse(J);
        mult(dadu, J, da);

        for (int i=0; i<N; i++){
            tet.a[i][j]   = a[i];
            tet.dadx[i][j]= da(i,0);
            tet.dady[i][j]= da(i,1);
            tet.dadz[i][j]= da(i,2);
            }
        tet.weight[j]    = detJ * pds[j];
    }
}

}
