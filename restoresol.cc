#include "fem.h"

#include "pt3D.h"

using namespace std;

void Fem::restoresol(double scaling, string *filename)  // filename may be NULL
{
string str("sol.in");

if (filename){
    str = *filename;
    }

ifstream fin(str, std::ifstream::in); //  *ct*
if (!fin){
    if(VERBOSE) { cerr << "pb ouverture fichier: " << str << "en lecture" << endl; }
    SYSTEM_ERROR;}

getline(fin, str); // 1eme ligne

unsigned long idx = str.find(":");
idx +=2;

t = stod(str.substr(idx));

if(VERBOSE) { cout << "fichier solution: " << str << " a l'instant t = " << t << endl; }

for (int i=0; i<NOD; i++){
    Node &n = node[i];
    Node node_;
    int i_;
    fin >> i_ >>  node_.p >> n.u >> n.phi;// carefull! >> is overloaded for class pt3D

    node_.p *= scaling;

    if (( Pt::norme2(n.p - node_.p) > sq(diam * 1e-9))&&VERBOSE) 
	{
    cerr << "WARNING difference dans la position des noeuds"<< endl;
    cerr << i  << "\t" << n.p << endl << i_ << "\t" << node_.p << endl;
    }
    
    if (i!=i_){
        if(VERBOSE) { cerr << "fichier .sol incompatibilite de noeuds"<< endl; }
        SYSTEM_ERROR;
        }
    }

fin.close();    			     
}
