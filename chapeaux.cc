#include "linear_algebra.h"

using namespace std;

void Fem::chapeaux(void)
{
//const int FAC = fem.FAC;
//const int TET = fem.TET;

/********************* FACES *******************/
for (int i_f=0; i_f<FAC; i_f++){
   Facette::Fac &f = fac[i_f];
   
//const int NPI = Fac::NPI;// NPI = 4

   //double u[Facette::NPI]   = {   1/3.,   1/5.,   3/5.,   1/5.};
   //double v[Facette::NPI]   = {   1/3.,   1/5.,   1/5.,   3/5.};
   //double pds[Facette::NPI] = {-27/96., 25/96., 25/96., 25/96.};
  
   double detJ = 2*f.surf;

   for (int j=0; j<Facette::NPI; j++){
       f.a[0][j] = 1.-Facette::u[j]-Facette::v[j];
       f.a[1][j] = Facette::u[j];
       f.a[2][j] = Facette::v[j];
       f.weight[j]    = detJ * Facette::pds[j];
       }
   }

/****************** TETRAS *****************/
for (int i_t=0; i_t<TET; i_t++){
    Tetra::Tet &t = tet[i_t];
    
gmm::dense_matrix <double> nod(3,Tetra::N); 
    for (int i=0; i<Tetra::N; i++){
        int i_= t.ind[i];
//	cout << "t:" << t <<",  nod:" << i_ << endl;
        nod(0,i) = node[i_].x;
        nod(1,i) = node[i_].y;
        nod(2,i) = node[i_].z;
        }
// cout << nod <<endl;

    for (int j=0; j<Tetra::NPI; j++){
        vector <double> a(Tetra::N);
        gmm::dense_matrix <double> da(Tetra::N,3), dadu(Tetra::N,3), J(3,3);
        a[0] = 1.-Tetra::u[j]-Tetra::v[j]-Tetra::w[j];
        a[1] = Tetra::u[j];
        a[2] = Tetra::v[j];
        a[3] = Tetra::w[j];

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
            what << "Singular jacobian in tetrahedron " << i_t;
            throw runtime_error(what.str());
#else
            cerr << "jacobienne singuliere ds le tetraedre " << i_t << endl;
            exit(1);
#endif
            }

        lu_inverse(J);
        mult(dadu, J, da);

        for (int i=0; i<Tetra::N; i++){
            t.a[i][j]   = a[i];
            t.dadx[i][j]= da(i,0);
            t.dady[i][j]= da(i,1);
            t.dadz[i][j]= da(i,2);
            }
        t.weight[j]    = detJ * Tetra::pds[j];
    }
}

}
