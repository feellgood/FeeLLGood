#include "fem.h"

void savecfg_vtk(Fem &fem, int nt, string *filename)  // filename may be NULL
{
string str;
ostringstream ostr; 
if (filename){
    str = *filename;
    }
else{
    ostr.str("");
    ostr << fem.simname  << boost::format("_%d_B%6f_iter%d.vtk") % fem.SEQ % fem.Bext % nt;
    str = ostr.str();
    }
IF_VERBOSE(fem) {
cout <<"\n -------------------" << endl;
cout << " " << str << endl; 
}

ofstream fout(str.c_str(), ios::out);
if (!fout){
    IF_VERBOSE(fem) cerr << "pb ouverture fichier : " << str << "en ecriture" << endl;
    SYSTEM_ERROR;}

ostr.str("");
ostr << fem.simname  << boost::format(" time : %+20.10e") % fem.t;
str = ostr.str();

const int NOD = fem.NOD;
const int TET = fem.TET;
const int N   = Tet::N;

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
    fout << boost::format("%+20.10e %+20.10e %+20.10e") % x % y % z << endl;
    }

fout << boost::format("CELLS %8d %8d") % TET % (5*TET) << endl;
for (int t=0; t<TET; t++){
    Tet &tet = fem.tet[t];
    fout << setw(8) << N;
    for (int i=0; i<N; i++)
        fout << setw(8) << tet.ind[i];
    fout << endl;
    }

fout << boost::format("CELL_TYPES %8d") % TET << endl;
for (int t=0; t<TET; t++)
    fout << setw(8) << 10 << endl;
   
fout << boost::format("POINT_DATA %8d") % NOD << endl;
fout << boost::format("SCALARS phi float %d") % 1 << endl;
fout << "LOOKUP_TABLE my_table" << endl;

for (int i=0; i<NOD; i++){
    Node &node = fem.node[i];
    double u1  = node.u[0];
    double u2  = node.u[1];
    double u3  = node.u[2];
    double phi = node.phi;
    fout << boost::format("%+20.10e") % phi << endl;
    }

fout << "VECTORS u float" << endl;

for (int i=0; i<NOD; i++){
    Node &node = fem.node[i];
    double u1  = node.u[0];
    double u2  = node.u[1];
    double u3  = node.u[2];
    double phi = node.phi;
    fout << boost::format("%+20.10e %+20.10e %+20.10e") % u1 % u2 % u3 << endl;
    }
}    
