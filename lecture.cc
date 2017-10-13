#include "fem.h"
#define DEBUG 0

void lecture(Fem &fem, double scale, Regions *regions)  // regions may be NULL
{
int NOD, ELEM, tags, reg, TYP;
double value;
string trash, symb;
pair <string,int> p;

string str = fem.pbname;
ifstream msh(str.c_str());   // ouverture du fichier probleme en lecture
if (!msh){
    IF_VERBOSE(fem) cerr << "Impossible d'ouvrir le fichier " << str << endl;
    SYSTEM_ERROR;}

while(msh >> symb){
    if (symb == "$Scale"){
        /* scale == 0.0 is a sentinel to mean read from the file; floating point
         * comparison is reliable in this case */
        if (scale == 0.0) msh >> scale;      // lecture de l'echelle
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
    cerr << "could not find $Nodes" << endl;
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
fem.scale = scale;

msh >> NOD;        // lecture des noeuds
fem.NOD = NOD; 

fem.node.resize(NOD);
for (int i=0; i<NOD; i++){
    double x,y,z;
    msh >> trash >> x >> y >> z;

    fem.node[i].x = x * scale;
    fem.node[i].y = y * scale;
    fem.node[i].z = z * scale;
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
            Fac fac;
            fac.reg = reg;
	    
            msh >> fac.ind[0] >> fac.ind[1] >> fac.ind[2];
	    for (int i=0; i<3; i++)
	        fac.ind[i]--;           // passage convention Matlab/msh a C++
		
            fem.fac.push_back(fac);
            if (regions) ++regions->surfaces[reg];
            break;
	    }
        case 4:{
            Tet tet;
            tet.reg = reg; 
	    
            msh >> tet.ind[0] >> tet.ind[1] >> tet.ind[2] >> tet.ind[3];
	    for (int i=0; i<4; i++)
	    	tet.ind[i]--;           // passage convention Matlab/msh a C++
		
            fem.tet.push_back(tet);
            if (regions) ++regions->volumes[reg];
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

fem.FAC = fem.fac.size();
fem.TET = fem.tet.size();
fem.SRC = fem.FAC * Fac::NPI + fem.TET * Tet::NPI;

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
//cout << symb << '\t' << reg << '\t' << value << endl;
    p = make_pair(symb,reg);
    fem.param[p] = value;
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
