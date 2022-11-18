#include "mesh.h"
#include "feellgoodSettings.h"
#include "pt3D.h"


void on_fail_msg_error(std::ifstream &f_in, const std::string strWhat)
{
if (f_in.fail())
    {
    std::cerr << strWhat << std::endl; 
    SYSTEM_ERROR;
    } 
}


void lookFor(const bool verbose, std::ifstream &f_in, const std::string strWhat)
{
std::string symb = "";
while (symb != strWhat) {f_in >> symb;}

if (verbose) on_fail_msg_error(f_in,"could not find beacon " + strWhat);
}


void mesh::readMesh(Settings const& mySets)
{
if(mySets.verbose) { std::cout << "Reading mesh file " << mySets.getPbName() << ":\n"; }
std::string symb;
std::ifstream msh( mySets.getPbName() );

on_fail_msg_error(msh,"cannot open file");

msh >> symb;
if(symb == "$MeshFormat")
    {
    std::string mshFormat = "";
    msh >> mshFormat;
    if(mshFormat == "2.2" ) 
        {
        if(mySets.verbose) {std::cout << "  file format: 2.2\n";}
        int tags, reg, TYP;
  
    lookFor(mySets.verbose,msh,"$Nodes");
    
    double scale = mySets.getScale();
    int nbNod;
    msh >> nbNod;
    init_node(nbNod);

    for (int i=0; i<nbNod; i++)
        {
        msh >> symb >> node[i].p;
        node[i].p.rescale(scale);
        }

    on_fail_msg_error(msh,"error while reading nodes" );
    
    lookFor(mySets.verbose,msh,"$Elements");

    int nbElem;
    msh >> nbElem; 
    if(mySets.verbose) {std::cout << "  element count: " << nbElem << '\n';}
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
                tet.push_back( Tetra::Tet(node,reg,mySets.findTetraRegionIdx(reg),i0,i1,i2,i3) );
                break;
                }
            default:
                std::cerr<< "unknown type in mesh $Elements" <<std::endl;
            break;
            }
        }

        on_fail_msg_error(msh,"error while reading elements" );
        }
    else { std::cout <<"mesh file format " << mshFormat << " not supported." << std::endl; SYSTEM_ERROR; }
    }

if (! msh.fail())
    {
    if(mySets.verbose) std::cout<< "  closing file\n";
    msh.close();
    }
else
    { std::cout<< "error before closing mesh." <<std::endl; SYSTEM_ERROR; }
}

double mesh::readSol(bool VERBOSE,double scaling, const std::string fileName)
{
std::ifstream fin(fileName, std::ifstream::in);

on_fail_msg_error(fin, "cannot open file: " + fileName );

std::string str;
getline(fin, str);

unsigned long idx = str.find(":");
idx +=2;

double t = stod(str.substr(idx));

if(VERBOSE) { std::cout << ".sol file: " << str << " @ time t = " << t << std::endl; }

for (unsigned int i=0; i<node.size(); i++)
    {
    Nodes::Node & n = node[i];
    Nodes::Node node_;
    unsigned int i_;
    fin >> i_ >>  node_.p >> n.u >> n.phi;// carefull! >> is overloaded for class pt3D

    if (i!=i_)
        {
        if(VERBOSE) { std::cerr << ".sol file mismatch with mesh nodes"<< std::endl; }
        fin.close();
        SYSTEM_ERROR;
        }
    else
        {
        node_.p.rescale(scaling);
        if (( Pt::norme2(n.p - node_.p) > Pt::sq(diam * scaling))&&VERBOSE) 
            { std::cerr << "WARNING difference in node positions"  << i  << "\t" << n.p << std::endl << i_ << "\t" << node_.p << std::endl; }
        }
    }
fin.close();    			     

return t;
}

