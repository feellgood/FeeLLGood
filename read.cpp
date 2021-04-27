#include "mesh.h"

#include "feellgoodSettings.h"

#include "pt3D.h"

using namespace std;


void mesh::readOldMesh(Settings const& mySets,ifstream &msh)
{
    int tags, reg, TYP;
    string symb = "";
   
while (symb != "$Nodes") {msh >> symb;}

if (msh.fail())
    {
    if(mySets.verbose) { cerr << "could not find $Nodes" << endl; }
    SYSTEM_ERROR;
    }

double scale = mySets.getScale();

msh >> nbNod;
init_node(nbNod);

for (int i=0; i<nbNod; i++)
    {
    msh >> symb >> node[i].p;
	node[i].p.rescale(scale);
	}

if (msh.fail())
    {
    if(mySets.verbose) {cerr << "error while reading nodes" << endl;}
    SYSTEM_ERROR;
    }

while (symb != "$Elements") {msh >> symb;}

if (msh.fail())
    {
    if(mySets.verbose) {std::cerr << "could not find $Elements" << std::endl;}
    SYSTEM_ERROR;
    }
    
int nbElem;
msh >> nbElem; 
if(mySets.verbose) {std::cout <<"nb Elem: "<< nbElem <<std::endl;}
while((msh >> symb)&&(symb != "$EndElements")&&(! msh.fail() ) ) 
    {
    msh >> TYP >> tags >> reg;
    for (int i=1; i<tags; i++)
        msh >> symb;
    switch (TYP){
        case 2:{
            int i0,i1,i2;
            msh >> i0 >> i1 >> i2;
            fac.push_back( Facette::Fac(node,nbNod,reg,mySets.findFacetteRegionIdx(reg),i0,i1,i2 ) );
            break;
	    }
        case 4:{
            int i0,i1,i2,i3;
            msh >> i0 >> i1 >> i2 >> i3;
            tet.push_back( Tetra::Tet(node,nbNod,reg,mySets.findTetraRegionIdx(reg),i0,i1,i2,i3) );
            break;
	    }
        default:
            std::cerr<< "unknown type in mesh $Elements" <<std::endl;
        break;
        }
    }

if(mySets.verbose) { std::cout << "last symb = " << symb << std::endl; }
if ((symb != "$EndElements") && msh.fail()) 
    {std::cerr << "error while reading elements; symb = " << symb << std::endl;SYSTEM_ERROR;}
}

void mesh::readNewMesh(Settings const& mySets,ifstream &msh)
{
string symb = "";
   
while (symb != "$Nodes") {msh >> symb;}

if (msh.fail())
    {
    if(mySets.verbose) { cerr << "could not find $Nodes" << endl; }
    SYSTEM_ERROR;
    }

double scale = mySets.getScale();

int numBloc,idx_begin,idx_end;
msh >> numBloc >> nbNod >> idx_begin >> idx_end;

if (idx_end != nbNod) {std::cout << "indexation not supported" << std::endl;SYSTEM_ERROR;}

init_node(nbNod);

int i=0;
while (i<nbNod)
    {
    int j=0;
    while(j<numBloc)
        {int TYP,nbNodBloc;
        
        msh >> TYP >> symb >> idx_begin >> idx_end;
        if (TYP == 0)
            {msh >>symb;msh >> node[i].p; node[i].p.rescale(scale);i++;}
        else if ((TYP == 1)||(TYP == 2)||(TYP == 3))
            {nbNodBloc = idx_end;
            std::cout<< "bloc#"<< j <<std::endl;
            for(int k=0;k<nbNodBloc;k++) {msh >> symb;std::cout << "symb=" << symb << std::endl;} // reading of the nodes index, we don't care
            for(int k=0;k<nbNodBloc;k++) {msh >> node[i].p;std::cout << "node=" << node[i].p << std::endl; node[i].p.rescale(scale);i++;} //x,y,z of the nodes
            }
        else {std::cout<<"oulala"<<std::endl;SYSTEM_ERROR;}
        j++;
        }
    }

    msh >> symb;
    std::cout << "symb=" << symb << std::endl;
    if (symb != "$EndNodes") {std::cout << "error while reading nodes." << std::endl;SYSTEM_ERROR;}

while (symb != "$Elements") {msh >> symb;}

if (msh.fail())
    {
    if(mySets.verbose) {cerr << "could not find $Elements" << endl;}
    SYSTEM_ERROR;
    }
    
    std::cout << "to be continued..." << std::endl;
    SYSTEM_ERROR;
int TYP,tags,reg;
while ((symb != "$EndElements")&&(symb != "$End")&& !(msh.fail()) )
    {
    msh >> symb; 
    msh >> TYP >> tags >> reg;
    for (int i=1; i<tags; i++)
        msh >> symb;
    switch (TYP){
        case 2:{
            int i0,i1,i2;
            msh >> i0 >> i1 >> i2;
            
            fac.push_back( Facette::Fac(node,nbNod,reg,mySets.findFacetteRegionIdx(reg),i0,i1,i2 ) );
            break;
	    }
        case 4:{
            int i0,i1,i2,i3;
            msh >> i0 >> i1 >> i2 >> i3;
            
            tet.push_back( Tetra::Tet(node,nbNod,reg,mySets.findTetraRegionIdx(reg),i0,i1,i2,i3) );
            break;
	    }
        default:
            std::cout<< "unknown type in mesh $Elements" <<std::endl;
        break;
        }
    }

if ((symb != "$EndElements") && msh.fail()) 
    {cerr << "error while reading elements; symb = " << symb << endl;SYSTEM_ERROR;}
}

void mesh::readMesh(Settings const& mySets)
{
string symb;
ifstream msh( mySets.getPbName() );

if (msh.fail())
    {
    cerr << "cannot open file " << mySets.getPbName() << endl; 
    SYSTEM_ERROR;
    } 

msh >> symb;
if(symb == "$MeshFormat")
    {
    string mshFormat = "";
    msh >> mshFormat;
    if(mshFormat == "2.2" ) 
        {
        if(mySets.verbose) {std::cout << "mesh file format 2.2" << std::endl;} 
        readOldMesh(mySets,msh);
        }
    else if (mshFormat == "4.0")
        {
        if(mySets.verbose)  {std::cout <<"mesh file format 4.0 not supported, only 2.2 and 4.1 are feellgood readable."<< std::endl; } 
        SYSTEM_ERROR;
        }
    else if (mshFormat == "4.1") 
        {
        if (mySets.verbose) {std::cout <<"mesh file format 4.1 not supported yet, coming soon" << std::endl; }
        SYSTEM_ERROR;
        readNewMesh(mySets,msh); 
        }
    else { std::cout <<"mesh file format " << mshFormat << " not supported." << std::endl; SYSTEM_ERROR; }
    }

if (! msh.fail())
    {
    if(mySets.verbose) std::cout<< "closing mesh" <<std::endl;
    msh.close();
    }
else
    { std::cout<< "error while closing mesh." <<std::endl; SYSTEM_ERROR; }
}

double mesh::readSol(bool VERBOSE,double scaling, const string fileName)
{
ifstream fin(fileName, std::ifstream::in);
if (!fin)
    {
    cerr << "cannot open file: " << fileName << endl;
    SYSTEM_ERROR;
    }

string str;
getline(fin, str);

unsigned long idx = str.find(":");
idx +=2;

double t = stod(str.substr(idx));

if(VERBOSE) { cout << ".sol file: " << str << " @ time t = " << t << endl; }

const int NOD = getNbNodes();

for (int i=0; i<NOD; i++)
    {
    Nodes::Node & n = node[i];
    Nodes::Node node_;
    int i_;
    fin >> i_ >>  node_.p >> n.u >> n.phi;// carefull! >> is overloaded for class pt3D

    if (i!=i_)
        {
        if(VERBOSE) { cerr << ".sol file mismatch with mesh nodes"<< endl; }
        fin.close();
        SYSTEM_ERROR;
        }
    else
        {
        node_.p.rescale(scaling);
        if (( Pt::norme2(n.p - node_.p) > Pt::sq(diam * scaling))&&VERBOSE) 
            { cerr << "WARNING difference in node positions"  << i  << "\t" << n.p << endl << i_ << "\t" << node_.p << endl; }
        }
    }
fin.close();    			     

return t;
}

