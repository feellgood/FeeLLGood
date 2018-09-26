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
if(VERBOSE) { cout << " " << str << endl; }

ofstream fout(str, ios::out);
if (!fout){
   if(VERBOSE) { cerr << "pb ouverture fichier " << str << "en ecriture" << endl; }
   SYSTEM_ERROR;}
//fout << boost::format("#time : %+20.10e ") % fem.t << endl;
fout << "#time : " << t <<endl;

for (int i=0; i<NOD; i++){
    Node &n = node[i];
    Pt::pt3D p = n.p/s;
    //   fout << boost::format("%8d %+20.10f %+20.10f %+20.10f %+20.10f %+20.10f %+20.10f %+20.10e") 
//                   % i % x % y % z % u1 % u2 % u3 % phi << endl;}
	fout << i << "\t" << p << "\t" << n.u << "\t" << n.phi << endl; //carefull! << overloaded
	}

fout.close();
}
