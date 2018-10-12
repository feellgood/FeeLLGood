#include <algorithm>

#include "fem.h"
#include "config.h"

//#include <boost/format.hpp>


using namespace std;


void Fem::saver(Settings &settings, ofstream &fout, int nt)
{
int n1 = settings.n1;
int n2 = settings.n2;

if ((nt%n1)==0)
    {
   // fout << boost::format("%8d %+20.10e %+20.10e %+20.10e") % nt % t % dt % dumax;
   // fout << boost::format("%+20.10e %+20.10e %+20.10e") % u_moy(fem,0) % u_moy(fem,1) % u_moy(fem,2);
   // fout << boost::format("%+20.10e %+20.10e %+20.10e %+20.10e %+20.10e %+20.10e %+20.10e")% Ee % Ea % Ed % Ez % Etot % DW_z % DW_vz<< endl;
    fout << nt <<"\t" << t <<"\t" << settings.dt <<"\t" << vmax*settings.dt <<"\t";
    fout << moy<U>(Pt::IDX_X) <<"\t" << moy<U>(Pt::IDX_Y) <<"\t" << moy<U>(Pt::IDX_Z) <<"\t";
    fout << E[0] <<"\t" << E[1] <<"\t" << E[2] <<"\t" << E[3] <<"\t" << Etot <<"\t" << DW_z <<"\t" << DW_vz <<"\t" << endl;
    }

    string baseName = settings.r_path_output_dir + settings.getSimName();
    
if ((nt%n2)==0) 
    {
    if (settings.withVtk)
        {
        string str = baseName + "_" + to_string(SEQ) + "_B" + to_string(Bext) + "_iter" + to_string(nt) + ".vtk";
            //boost::format("_%d_B%6f_iter%d.vtk") % fem.SEQ % fem.Bext % nt;
        savecfg_vtk(str);
        }
    
    string str = baseName + "_" + to_string(SEQ) +"_B" + to_string(Bext) + "_iter" + to_string(nt) + ".sol";
 //<< boost::format("_%d_B%6f_iter%d.sol") % fem.SEQ % fem.Bext % nt;
    savesol(str,settings.getScale());
    
    //string str = baseName + "_" + to_string(SEQ) + "_B" + to_string(Bext) + "_iter" + to_string(nt) + ".hdm";
    //    saveH(str);
    }
}

void Fem::savecfg_vtk(string fileName)
{
    const int TET = tet.size();
    
if(VERBOSE) { cout <<"\n -------------------\n " << fileName << endl; }

ofstream fout(fileName, ios::out);
if (!fout){
    if(VERBOSE) cerr << "cannot open file : " << fileName << endl;
    SYSTEM_ERROR;}


fout << "# vtk DataFile Version 2.0" << endl;
fout << "time : " << t << endl; // boost::format(" time : %+20.10e") % fem.t;
fout << "ASCII" << endl;
fout << "DATASET UNSTRUCTURED_GRID" << endl;
fout << "POINTS "<< NOD << " float" << endl;

for (int i=0; i<NOD; i++){
fout << node[i].p.x() << "\t" << node[i].p.y() << "\t" << node[i].p.z() + DW_z << endl;
//fout << boost::format("%+20.10e %+20.10e %+20.10e") % x % y % z << endl;
    }

fout << "CELLS " << setw(8) << TET << "\t" << setw(8) << (5*TET) << endl;

for (int t=0; t<TET; t++){
    fout << setw(8) << Tetra::N;
    for (int i=0; i<Tetra::N; i++) { fout << setw(8) << tet[t].ind[i]; }
    fout << endl;
    }

fout << "CELL_TYPES " << setw(8) << TET << endl;

for (int t=0; t<TET; t++) { fout << setw(8) << 10 << endl; }

fout <<"POINT_DATA " << setw(8) << NOD << endl;
fout <<"SCALARS phi float 1" << endl;
fout << "LOOKUP_TABLE my_table" << endl;

for (int i=0; i<NOD; i++) { fout << node[i].phi << endl; }

fout << "VECTORS u float" << endl;

for (int i=0; i<NOD; i++){
    //fout << boost::format("%+20.10e %+20.10e %+20.10e") % u1 % u2 % u3 << endl;
fout << node[i].u << endl;
    }
}

void Fem::savesol(string fileName,double s)
{
if(VERBOSE) { cout << " " << fileName << endl; }

ofstream fout(fileName, ios::out);
if (!fout){
   if(VERBOSE) { cerr << "cannot open file " << fileName << endl; }
   SYSTEM_ERROR;}
//fout << boost::format("#time : %+20.10e ") % fem.t << endl;
fout << "#time : " << t <<endl;

//   fout << boost::format("%8d %+20.10f %+20.10f %+20.10f %+20.10f %+20.10f %+20.10f %+20.10e") 
//                   % i % x % y % z % u1 % u2 % u3 % phi << endl;}

int i=0;
std::for_each(node.begin(),node.end(),
    [&i,s,&fout](Node &n)
        { 
        Pt::pt3D p = n.p/s;
        fout << i << "\t" << p << "\t" << n.u << "\t" << n.phi << endl;
        i++;
        } // lambda works with << overloaded
    );
if(VERBOSE) { cout << i <<" nodes written." << endl; }

fout.close();
}

void Fem::saveH(string fileName,double scale)
{
const int TET = tet.size();

cout << " " << fileName << endl <<" -------------------" << endl << endl;

ofstream fout(fileName, ios::out);
if (!fout){
   cerr << "cannot open file " << fileName << endl;
   SYSTEM_ERROR;}
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

    tiny::mult<double, 3, Tetra::N, Tetra::NPI> (nod, Tetra::a, gauss);
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
