#include "linear_algebra.h"

using namespace std;

void Fem::chapeaux(double epsilon)
{
/********************* FACES *******************/
std::for_each(fac.begin(),fac.end(),[](Facette::Fac &f) {f.init();});


/****************** TETRAS *****************/
gmm::dense_matrix <double> dadu(Tetra::N,3);
Tetra::init_dadu(dadu);

std::for_each(tet.begin(),tet.end(),
    [this,epsilon,&dadu](Tetra::Tet &t)
    {
    gmm::dense_matrix <double> nod(3,Tetra::N); 
    t.getNod(nod,node);
    gmm::dense_matrix <double> J(3,3),da(Tetra::N,3);
    mult(nod, dadu, J);
    double detJ = lu_det(J);
    if (fabs(detJ) < epsilon){
        #ifdef LIBRARY
            ostringstream what;
            what << "Singular jacobian in tetrahedron ";
            throw runtime_error(what.str());
        #else
            cerr << "jacobienne singuliere ds le tetraedre " << endl;
			t.infos();
            exit(1);
        #endif
            }
    lu_inverse(J);
    mult(dadu, J, da);
    
    for (int j=0; j<Tetra::NPI; j++)
        {
        for (int i=0; i<Tetra::N; i++)
            { t.dadx[i][j]=da(i,0); t.dady[i][j]=da(i,1); t.dadz[i][j]=da(i,2); }
        t.weight[j]    = detJ * Tetra::pds[j];
        }    
    }
);//end for_each sur tet 
}
