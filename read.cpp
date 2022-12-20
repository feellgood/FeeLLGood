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


bool lookFor(const bool _b, std::ifstream &f_in, const std::string strWhat)
{
std::string symb = "";
while ( (f_in.peek()!=EOF) && (symb != strWhat)) {f_in >> symb;}

if (_b) on_fail_msg_error(f_in,"could not find beacon " + strWhat);

return !(f_in.fail());
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


	bool beaconFound = lookFor(false,msh,"$PhysicalNames");
	
	int nbRegNames;
	if (beaconFound)
		{
		msh >> nbRegNames;
		if(mySets.verbose) { std::cout << "found $PhysicalNames beacon : " << nbRegNames <<" reg names.\n"; }
		
		while((msh >> symb)&&(symb != "$EndPhysicalNames")&&(! msh.fail() ) ) 
        		{
        		std::string name;
        		msh >> tags  >> name;
        		
        		switch( stoi(symb) ){
        			case 2:{
        				surfRegNames[name] = tags;
        			break;}
        			case 3:{
        				volRegNames[name] = tags;
        			break;}
        			default:
        				std::cerr<< "unknown type in mesh $PhysicalNames" <<std::endl;
        			break;
        			}
        		}
        	if(mySets.verbose)
        		{
        		std::map<std::string, int>::iterator it;
   			for(it=surfRegNames.begin(); it!=surfRegNames.end(); ++it){ std::cout << it->first << " => " << it->second << '\n';}
        		for(it=volRegNames.begin(); it!=volRegNames.end(); ++it){ std::cout << it->first << " => " << it->second << '\n';}
        		}
        	
		}
  	else
  		{
  		if (mySets.verbose)
  			{ std::cout << "No $PhysicalNames defined.\n"; }
  		
  		msh.clear(); 
  		msh.seekg(std::ios::beg); 
  		}
  
    beaconFound = lookFor(true,msh,"$Nodes");
    
    double scale = mySets.getScale();
    int nbNod;
    
    if (beaconFound)
        {
        msh >> nbNod;
        init_node(nbNod);

        for (int i=0; i<nbNod; i++)
            {
            msh >> symb >> node[i].p;
            node[i].p.rescale(scale);
            }
        }
        
    on_fail_msg_error(msh,"error while reading nodes" );
    
    beaconFound = lookFor(true,msh,"$Elements");

    int nbElem;
    
    if (beaconFound)
    {
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

