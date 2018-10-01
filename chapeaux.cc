#include "linear_algebra.h"

void Fem::chapeaux(double epsilon)
{
/********************* FACES *******************/
std::for_each(fac.begin(),fac.end(),[](Facette::Fac &f) {f.init();});


/****************** TETRAS *****************/
std::for_each(tet.begin(),tet.end(),
    [this,epsilon](Tetra::Tet &t)
    {
    double J[Pt::DIM][Pt::DIM];
    double detJ = t.Jacobian(J,node);
    double da[Tetra::N][Pt::DIM];
    
    if (fabs(detJ) < epsilon){
        #ifdef LIBRARY
            ostringstream what;
            what << "Singular jacobian in tetrahedron ";
            throw runtime_error(what.str());
        #else
            std::cerr << "jacobienne singuliere ds le tetraedre " << std::endl;
			t.infos();
            SYSTEM_ERROR;
        #endif
            }
    Pt::inverse(J,detJ);
    tiny::mult<double, Tetra::N, Pt::DIM, Pt::DIM> (Tetra::dadu, J, da);
    
    for (int j=0; j<Tetra::NPI; j++)
        {
        for (int i=0; i<Tetra::N; i++)
            { t.dadx[i][j]=da[i][0]; t.dady[i][j]=da[i][1]; t.dadz[i][j]=da[i][2]; }
        t.weight[j]    = detJ * Tetra::pds[j];
        }    
    }
);//end for_each sur tet 
}
