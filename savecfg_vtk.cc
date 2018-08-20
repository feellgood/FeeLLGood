#include "fem.h"

using namespace std;

void savecfg_vtk(Fem &fem,string baseName,double s, int nt, string *filename)  // filename may be NULL
{
string str;

if (filename) { str = *filename; }
else{
    str = baseName + "_" + to_string(fem.SEQ) + "_B" + to_string(fem.Bext) + "_iter" + to_string(nt) + ".vtk";
 //<< boost::format("_%d_B%6f_iter%d.vtk") % fem.SEQ % fem.Bext % nt;
    }

IF_VERBOSE(fem) { cout <<"\n -------------------" << endl << " " << str << endl; }

ofstream fout(str, ios::out);
if (!fout){
    IF_VERBOSE(fem) cerr << "pb ouverture fichier : " << str << "en ecriture" << endl;
    SYSTEM_ERROR;}

str = baseName + " time : " + to_string(fem.t);  
// << boost::format(" time : %+20.10e") % fem.t;

const int NOD = fem.NOD;
const int TET = fem.TET;

fout << "# vtk DataFile Version 2.0" << endl;
fout << str << endl;
fout << "ASCII" << endl;
fout << "DATASET UNSTRUCTURED_GRID" << endl;
fout << "POINTS "<< NOD << " float" << endl;

for (int i=0; i<NOD; i++){
    Node &node = fem.node[i];
    double x   = node.x;
    double y   = node.y;
    double z   = node.z+fem.DW_z;
fout << x << "\t" << y << "\t" << z << endl;    //fout << boost::format("%+20.10e %+20.10e %+20.10e") % x % y % z << endl;
    }

fout << "CELLS " << setw(8) << TET << "\t" << setw(8) << (5*TET) << endl;

for (int t=0; t<TET; t++){
    Tetra::Tet &tet = fem.tet[t];
    fout << setw(8) << Tetra::N;
    for (int i=0; i<Tetra::N; i++)
        fout << setw(8) << tet.ind[i];
    fout << endl;
    }

fout << "CELL_TYPES " << setw(8) << TET << endl;

for (int t=0; t<TET; t++) { fout << setw(8) << 10 << endl; }

fout <<"POINT_DATA " << setw(8) << NOD << endl;
fout <<"SCALARS phi float 1" << endl;
   
//fout << boost::format("POINT_DATA %8d") % NOD << endl;
//fout << boost::format("SCALARS phi float %d") % 1 << endl;

fout << "LOOKUP_TABLE my_table" << endl;

for (int i=0; i<NOD; i++){
    Node &node = fem.node[i];
    double phi = node.phi;
    //fout << boost::format("%+20.10e") % phi << endl;
fout << phi << endl;    
}

fout << "VECTORS u float" << endl;

for (int i=0; i<NOD; i++){
    Node &node = fem.node[i];
    double u1  = node.u[0];
    double u2  = node.u[1];
    double u3  = node.u[2];
    //fout << boost::format("%+20.10e %+20.10e %+20.10e") % u1 % u2 % u3 << endl;
fout << u1 << "\t" << u2 << "\t" << u3 << endl;
    }
}    
