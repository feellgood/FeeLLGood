#include "fem.h"

using namespace std;

void Fem::savesol(string baseName,double s, int nt, string *filename)
{
string str;

if (filename) { str = *filename; }
else{
    str = baseName + "_" + to_string(SEQ) +"_B" + to_string(Bext) + "_iter" + to_string(nt) + ".sol";
 //<< boost::format("_%d_B%6f_iter%d.sol") % fem.SEQ % fem.Bext % nt;
    }
IF_VERBOSE() cout << " " << str << endl;

ofstream fout(str, ios::out);
if (!fout){
   IF_VERBOSE() cerr << "pb ouverture fichier " << str << "en ecriture" << endl;
   SYSTEM_ERROR;}
//fout << boost::format("#time : %+20.10e ") % fem.t << endl;
fout << "#time : " << t <<endl;

for (int i=0; i<NOD; i++){
    Node &n = node[i];
    double x   = n.p.x() / s;
    double y   = n.p.y() / s;
    double z   = n.p.z() / s;
    double u1  = n.u[0];
    double u2  = n.u[1];
    double u3  = n.u[2];
    double phi = n.phi;
 
//   fout << boost::format("%8d %+20.10f %+20.10f %+20.10f %+20.10f %+20.10f %+20.10f %+20.10e") 
//                   % i % x % y % z % u1 % u2 % u3 % phi << endl;}
	fout << i << "\t" << x << "\t" << y << "\t" << z << "\t";
	fout << u1 << "\t" << u2 << "\t" << u3 << "\t" << phi << endl;
	}

fout.close();
}
