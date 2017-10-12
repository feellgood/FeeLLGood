#include "fem.h"
#include "tiny.h"

void saveH(Fem &fem, int nt)
{
string str;
str = fem.simname + "_" + to_string(fem.SEQ) + "_B" + to_string(fem.Bext) + "_iter" + to_string(nt) + ".hdm";

cout << " " << str << endl <<" -------------------" << endl << endl;

ofstream fout(str, ios::out);
if (!fout){
   cerr << "pb ouverture fichier " << str << "en ecriture" << endl;
   exit(1);}
fout << "#time : "<< fem.t << endl;

const int    TET   = fem.TET;
const double scale = fem.scale;

for (int t=0; t<TET; t++){
    Tet &tet = fem.tet[t];
    const int N   = Tet::N;
    const int NPI = Tet::NPI;
    
    /*------------------- INTERPOL --------------------*/
    double nod[3][N], gauss[3][NPI];
    double negphi_nod[N], Hdx[NPI], Hdy[NPI], Hdz[NPI];

    for (int i=0; i<N; i++){
        int i_= tet.ind[i];
	    Node &node = fem.node[i_];
	    nod[0][i] = node.x;
	    nod[1][i] = node.y;
	    nod[2][i] = node.z;
		negphi_nod[i] = -node.phi;
	    }

    tiny::mult<double, 3, N, NPI> (nod, tet.a, gauss);
	tiny::transposed_mult<double, N, NPI> (negphi_nod, tet.dadx, Hdx);
	tiny::transposed_mult<double, N, NPI> (negphi_nod, tet.dady, Hdy);
	tiny::transposed_mult<double, N, NPI> (negphi_nod, tet.dadz, Hdz);
    /*---------------------------------------------------*/
    
    for (int npi=0; npi<NPI; npi++){
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
