#include "fem.h"

#include "pt3D.h"

using namespace std;


void Fem::readOldMesh(Settings const& mySets,ifstream &msh)
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
int NOD;    
msh >> NOD;
node.resize(NOD);

for (int i=0; i<NOD; i++)
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
    if(mySets.verbose) {cerr << "could not find $Elements" << endl;}
    SYSTEM_ERROR;
    }

msh >> symb; 

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
            
            fac.push_back( Facette::Fac(NOD,reg,mySets.findFacetteRegionIdx(reg),i0,i1,i2 ) );
            break;
	    }
        case 4:{
            int i0,i1,i2,i3;
            msh >> i0 >> i1 >> i2 >> i3;
            
            tet.push_back( Tetra::Tet(NOD,reg,mySets.findTetraRegionIdx(reg),i0,i1,i2,i3) );
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

void Fem::readNewMesh(Settings const& mySets,ifstream &msh)
{
    std::cout <<"mesh file format 4.1 reading function coming soon."<< std::endl;
    SYSTEM_ERROR;
}

void Fem::readMesh(Settings const& mySets)
{
string symb;
ifstream msh( mySets.getPbName() );

if (!msh)
    {
    if(mySets.verbose) { cerr << "cannot open file " << mySets.getPbName() << endl; } 
    SYSTEM_ERROR;
    }

msh >> symb;
if(symb == "$MeshFormat")
    {
    string mshFormat;
    msh >> mshFormat;
    if(mshFormat == "2.2" ) { std::cout << "mesh file format 2.2" << std::endl; readOldMesh(mySets,msh); }
    else if (mshFormat == "4.0") { std::cout <<"mesh file format 4.0 not supported, only 2.2 and 4.1 are feellgood readable."<< std::endl; SYSTEM_ERROR;}
    else if (mshFormat == "4.1") { std::cout <<"mesh file format 4.1" << std::endl; readNewMesh(mySets,msh); }
    }
        
msh.close();
}

void Fem::readSol(bool VERBOSE,double scaling, string fileName)
{
ifstream fin(fileName, std::ifstream::in);
if (!fin){
    if(VERBOSE) { cerr << "cannot open .sol file: " << fileName << endl; }
    SYSTEM_ERROR;}

string str;
getline(fin, str); // 1eme ligne

unsigned long idx = str.find(":");
idx +=2;

t = stod(str.substr(idx));

if(VERBOSE) { cout << ".sol file: " << str << " @ time t = " << t << endl; }

const int NOD = node.size();

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
        if (( Pt::norme2(n.p - node_.p) > sq(diam * 1e-9))&&VERBOSE) 
            { cerr << "WARNING difference in node positions"  << i  << "\t" << n.p << endl << i_ << "\t" << node_.p << endl; }
        }
    }
fin.close();    			     
}

