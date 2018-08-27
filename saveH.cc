#include "fem.h"

using namespace std;

void Fem::saveH(string baseName,double scale, int nt)
{
string str = baseName + "_" + to_string(SEQ) + "_B" + to_string(Bext) + "_iter" + to_string(nt) + ".hdm";

cout << " " << str << endl <<" -------------------" << endl << endl;

ofstream fout(str, ios::out);
if (!fout){
   cerr << "pb ouverture fichier " << str << "en ecriture" << endl;
   exit(1);}
fout << "#time : "<< t << endl;

for (int t=0; t<TET; t++){
    Tetra::Tet &te = tet[t];
    
    /*------------------- INTERPOL --------------------*/
    double nod[3][Tetra::N], gauss[3][Tetra::NPI];
    double negphi_nod[Tetra::N], Hdx[Tetra::NPI], Hdy[Tetra::NPI], Hdz[Tetra::NPI];

    for (int i=0; i<Tetra::N; i++){
        int i_= te.ind[i];
	    Node &n = node[i_];
	    nod[0][i] = n.p.x();
	    nod[1][i] = n.p.y();
	    nod[2][i] = n.p.z();
		negphi_nod[i] = -n.phi;
	    }

    tiny::mult<double, 3, Tetra::N, Tetra::NPI> (nod, te.a, gauss);
	tiny::transposed_mult<double, Tetra::N, Tetra::NPI> (negphi_nod, te.dadx, Hdx);
	tiny::transposed_mult<double, Tetra::N, Tetra::NPI> (negphi_nod, te.dady, Hdy);
	tiny::transposed_mult<double, Tetra::N, Tetra::NPI> (negphi_nod, te.dadz, Hdz);
    /*---------------------------------------------------*/
    
    for (int npi=0; npi<Tetra::NPI; npi++){
	    double x = gauss[0][npi]/scale;
	    double y = gauss[1][npi]/scale;
	    double z = gauss[2][npi]/scale;
	
	
	fout << t << " " << npi << " " << x << " " << y << " " << z;
	fout << " "<< Hdx[npi] << " " << Hdy[npi] << " " << Hdz[npi] << endl;
	// " %8d %3d %+20.10f %+20.10f %+20.10f %+20.10e %+20.10e %+20.10e"
	}
}

fout.close();
}
