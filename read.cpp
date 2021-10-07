#include "mesh.h"

#include "feellgoodSettings.h"

#include "pt3D.h"

using namespace std;

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

