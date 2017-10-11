#include "fem.h"

void savesol(Fem &fem, int nt, string *filename)
{
string str;
ostringstream ostr; 
if (filename){
    str = *filename;
    }
else{
    ostr.str("");
    ostr << fem.simname  << boost::format("_%d_B%6f_iter%d.sol") % fem.SEQ % fem.Bext % nt;
    str = ostr.str();
    }
IF_VERBOSE(fem) cout << " " << str << endl;

ofstream fout(str.c_str(), ios::out);
if (!fout){
   IF_VERBOSE(fem) cerr << "pb ouverture fichier " << str << "en ecriture" << endl;
   SYSTEM_ERROR;}
fout << boost::format("#time : %+20.10e ") % fem.t << endl;

const int    NOD   = fem.NOD;
const double scale = fem.scale;

for (int i=0; i<NOD; i++){
    Node &node = fem.node[i];
    double x   = node.x / scale;
    double y   = node.y / scale;
    double z   = node.z / scale;
    double u1  = node.u[0];
    double u2  = node.u[1];
    double u3  = node.u[2];
    double phi = node.phi;
 
    fout << boost::format
    ("%8d %+20.10f %+20.10f %+20.10f %+20.10f %+20.10f %+20.10f %+20.10e") 
                   % i % x % y % z % u1 % u2 % u3 % phi << endl;}

fout.close();
}
