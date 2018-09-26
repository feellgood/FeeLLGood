#include "fem.h"

#include "pt3D.h"

#define DEBUG 0

using namespace std;

void Fem::lecture(Settings &mySets, double scale)//, Regions *regions)  // regions may be NULL
{
int ELEM, tags, reg, TYP;
double value;
string trash, symb;

string str = mySets.getPbName();
ifstream msh(str);   // ouverture du fichier probleme en lecture
if (!msh){
    if(VERBOSE) { cerr << "Impossible d'ouvrir le fichier " << str << endl; } 
    SYSTEM_ERROR;}

while(msh >> symb){
    if (symb == "$Scale"){
        /* scale == 0.0 is a sentinel to mean read from the file; floating point
         * comparison is reliable in this case */
        if (scale == 0.0) 
            {
            msh >> scale;      // lecture de l'echelle
            mySets.setScale( scale );
            }
        else msh >> trash;
        continue;
        }
    if (symb == "$Nodes")
        break;
    }

if (msh.fail()){
#ifdef LIBRARY
    throw runtime_error("could not find $Nodes");
#else
    if(VERBOSE) { cerr << "could not find $Nodes" << endl; }
    exit(1);
#endif
    }

if (scale == 0.0){
#ifdef LIBRARY
    throw runtime_error("scale must be passed as argument or defined in the problem file");
#else
    cerr << "scale must be passed as argument or defined in the problem file" << endl;
    exit(1);
#endif
    }

msh >> NOD;        // lecture des noeuds
node.resize(NOD);
for (int i=0; i<NOD; i++){
    msh >> trash >> node[i].p;
	node[i].p *= scale;
	//cout <<"scale=" << scale <<";" << trash << ";" << node[i].p << endl;
	}

if (msh.fail()){
#ifdef LIBRARY
    throw runtime_error("error while reading nodes");
#else
    cerr << "error while reading nodes" << endl;
    exit(1);
#endif
    }

while(msh >> symb){
    if (symb == "$Elements")
        break;
    }

if (msh.fail()){
#ifdef LIBRARY
    throw runtime_error("could not find $Elements");
#else
    cerr << "could not find $Elements" << endl;
    exit(1);
#endif
    }

msh >> ELEM;        // lecture des elements

while(msh >> symb){
    if (symb == "$EndElements" || (symb=="$End"))
        break;
    msh >> TYP >> tags >> reg;
    for (int i=1; i<tags; i++)
        msh >> trash;
    switch (TYP){
        case 2:{
            Facette::Fac f;
            f.reg = reg;
			f.idxPrm = mySets.findFacetteRegionIdx(reg);
            msh >> f.ind[0] >> f.ind[1] >> f.ind[2];
	    for (int i=0; i<3; i++)
	        f.ind[i]--;           // passage convention Matlab/msh a C++
		
            fac.push_back(f);
            //if (regions) ++regions->surfaces[reg];
            break;
	    }
        case 4:{
            Tetra::Tet te;
            te.reg = reg; 
			te.idxPrm = mySets.findTetraRegionIdx(reg);
            msh >> te.ind[0] >> te.ind[1] >> te.ind[2] >> te.ind[3];
	    for (int i=0; i<4; i++)
	    	te.ind[i]--;           // passage convention Matlab/msh a C++
		
            tet.push_back(te);
            //if (regions) ++regions->volumes[reg];
            break;
	    }
        default:
            getline(msh,trash);
	}
    }

if (msh.fail()){
#ifdef LIBRARY
    throw runtime_error("error while reading elements");
#else
    cerr << "error while reading elements" << endl;
    exit(1);
#endif
    }

FAC = fac.size();
TET = tet.size();
SRC = FAC * Facette::NPI + TET * Tetra::NPI;

while(msh >> symb){
    if (symb == "$Parameters")
        break;
    }

while(msh >> symb){  // lecture des parametres
    if (symb=="$EndParameters" || (symb=="$End"))
        break;
    if (msh.eof())  // be lenient
        break;
    msh >> reg >> value;
    if(VERBOSE) {cout <<"from msh: " << symb << '\t' << reg << '\t' << value << endl;}
    //p = make_pair(symb,reg);
    //mySets.param[p] = value;
    if (msh.fail()){
    #ifdef LIBRARY
        throw runtime_error("error while reading parameters");
    #else
        cerr << "error while reading parameters" << endl;
        exit(1);
    #endif
        }
    }

msh.close();
}
