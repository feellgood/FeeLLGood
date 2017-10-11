#include "fem.h"

void affichage(Fem &fem)
{
//cout << "process: " << getpid() <<"\n\n";

pair <string,int> p;
//map  <pair<string,int>,double> &param = fem.param; // unused *ct*

//p= make_pair("dt",-1);
//cout << "dt : " << param[p] << endl;

//p= make_pair("theta",-1);
//cout << "theta : " << param[p] << endl;

//p= make_pair("alpha",-1);
//cout << setw(8) << "alpha : " << param[p] << endl; 

int REG = fem.REG;
cout << "\n\t regions\t\t" << REG << endl;

int NOD = fem.NOD;
cout << "\t noeuds\t\t\t" << NOD << endl;
/*for (int np=0; np<NP; np++){
    Node &node=fem.node[np];
    cout<< setw(8) << node.ind << setw(8) << node.x << setw(8) << node.y << endl;      
    } */

int FAC = fem.FAC; 
cout << "\t faces\t\t\t" << FAC << endl;

int TET = fem.TET;
cout << "\t tetraedres\t\t" << TET << endl;
/*for (int ne=0; ne<NE; ne++){
    Element &e=fem.elt[ne];
    int TYP=e.TYP;
    int NBN=e.NBN;
    int NRG=e.NRG;
   cout<< setw(8) << TYP << setw(8) << NBN << setw(8) << NRG;
    for (int nbn=0; nbn<NBN; nbn++) {
        cout<< setw(8) << e.ind[nbn];
        }
    cout << endl;	
    } */

cout << "\t surface\t\t"  << fem.surf << endl;
cout << "\t volume\t\t\t" << fem.vol << endl;

/*
cout << "PROPRIETES_PHYSIQUES" << endl;

for (map<pair<string, int>, double>::const_iterator it=param.begin(); it!=param.end(); ++it)
  {
   p=it->first;
   cout<< p.first<< setw(8)<< p.second<< setw(8)<< it->second<< endl;
  }
*/
}
