#include <fstream>
#include <algorithm>
#include <gmsh.h>

#include "feellgoodSettings.h"
#include "mesh.h"
#include "meshUtils.h"

namespace Mesh
    {
void mesh::fromFile(Settings const &mySets)
    {
    gmsh::initialize();
    if (!mySets.verbose)
        { gmsh::option::setNumber("General.Terminal",0); } // to silent gmsh

    checkMeshFile(mySets);
    readNodes(mySets);
    readTetraedrons(mySets);
    for (unsigned int i = 0; i < tet.size(); i++)
        { tet[i].idx = i; }
    readTriangles(mySets);

    gmsh::clear();
    gmsh::finalize();
    }

void mesh::checkMeshFile(Settings const &mySets)
    {
    using namespace tags::msh;
    gmsh::open(mySets.getPbName());
    
    bool msh_available = (gmsh::model::getDimension() == DIM_OBJ_3D);
    if (!msh_available)
        { std::cout<<"Fatal Error: mesh file " << mySets.getPbName() << " is not 3D\n"; exit(1); }
    std::vector<int> elemTypes;
    gmsh::model::mesh::getElementTypes(elemTypes, DIM_OBJ_3D, -1);
    msh_available = msh_available && (elemTypes.size() == 1); // only one type of 3D element allowed
    msh_available = msh_available && (elemTypes[0] == TYP_ELEM_TETRAEDRON); // only Tetrahedrons
    if (!msh_available)
        { std::cout<<"Fatal Error: mesh file contains other 3D elements than tetraedrons.\n"; exit(1); }
    gmsh::model::mesh::getElementTypes(elemTypes, DIM_OBJ_2D, -1);
    msh_available = msh_available && (elemTypes.size() == 1); // only one type of 2D element allowed
    msh_available = msh_available && (elemTypes[0] == TYP_ELEM_TRIANGLE); // only Triangles
    if (!msh_available)
        { std::cout<<"Fatal Error: mesh file contains other 2D elements than triangles.\n"; exit(1); }
    bool o_2D = checkNamedObjects<Facette::prm>(mySets.paramFacette,DIM_OBJ_2D);
    bool o_3D = checkNamedObjects<Tetra::prm>(mySets.paramTetra,DIM_OBJ_3D);
    if (!(o_2D && o_3D))
        { std::cout <<"Fatal Error: mismatch between input mesh file and yaml regions.\n"; exit(1); }
    }

void mesh::readNodes(Settings const &mySets /**< [in] */)
    {
    std::vector<std::size_t> nodeT;
    std::vector<double> nodeC, nodeP;
    gmsh::model::mesh::getNodes(nodeT, nodeC, nodeP, -1,-1); // all nodes of the mesh
    node.resize( nodeT.size() );
    for(unsigned int i=0;i<nodeT.size();i++)
        {
        double scale = mySets.getScale();
        node[i].p = Eigen::Vector3d(nodeC[3*i],nodeC[3*i + 1],nodeC[3*i + 2]) * scale;
        }
    }

void mesh::readTetraedrons(Settings const &mySets /**< [in] */)
    {
    using namespace tags::msh;

    std::vector<std::pair<int, int> > entities;
    gmsh::model::getEntities(entities,DIM_OBJ_3D);
    std::for_each(entities.begin(),entities.end(), [this,&mySets]( std::pair<int, int> &p)
        {
        std::string entity_name;
        gmsh::model::getEntityName(p.first,p.second,entity_name);// warning! entity index p.second is different from physical region volume index
        std::vector<int> physicalTags;
        gmsh::model::getPhysicalGroupsForEntity(p.first, p.second, physicalTags);
        if(physicalTags.size())
            {
            for(auto physTag : physicalTags)
                {
                std::string n;
                gmsh::model::getPhysicalName(p.first, physTag, n);

                if ( std::any_of(mySets.paramTetra.begin(),mySets.paramTetra.end(),
                     [&n] (Tetra::prm const &p) { return p.regName == n; } ) )
                    { //named volume region is found
                    std::vector<int> elemTypes;
                    std::vector<std::vector<std::size_t> > elemTags, elemNodeTags;
                    gmsh::model::mesh::getElements(elemTypes, elemTags, elemNodeTags, p.first, p.second);

                    if (std::all_of(elemTypes.begin(), elemTypes.end(), [](int &e_type){return e_type == TYP_ELEM_TETRAEDRON;} ))
                        { // all elements are tetrahedrons
                        int idx = mySets.findTetraRegionIdx(n);
                        for(unsigned int i=0;i<elemNodeTags[0].size();i+=4)
                            {
                            int i0 = elemNodeTags[0][i];
                            int i1 = elemNodeTags[0][i+1];
                            int i2 = elemNodeTags[0][i+2];
                            int i3 = elemNodeTags[0][i+3];
                            if (idx > -1) tet.push_back(Tetra::Tet(node, idx, {i0,i1,i2,i3}));
                            }
                        }
                    else
                        { std::cout << "Fatal error: a volume mesh entity contains other elements than tetraedrons.\n"; exit(1); }
                    }
                }
            }
        } );// end for_each on entities
    }

void mesh::readTriangles(Settings const &mySets /**< [in] */)
    {
    using namespace tags::msh;

    std::vector<std::pair<int, int> > entities;
    gmsh::model::getEntities(entities,DIM_OBJ_2D);
    std::for_each(entities.begin(),entities.end(), [this,&mySets]( std::pair<int, int> &p)
        {
        std::string entity_name;
        gmsh::model::getEntityName(p.first,p.second,entity_name);// warning! entity index p.second is different from physical region surface index
        std::vector<int> physicalTags;
        gmsh::model::getPhysicalGroupsForEntity(p.first, p.second, physicalTags);
        if(physicalTags.size())
            {
            for(auto physTag : physicalTags)
                {
                std::string n;
                gmsh::model::getPhysicalName(p.first, physTag, n);
                if ( std::any_of(mySets.paramFacette.begin(),mySets.paramFacette.end(),
                     [&n] (Facette::prm const &p) { return p.regName == n; } ) )
                    { //named surface region is found
                    std::vector<int> elemTypes;
                    std::vector<std::vector<std::size_t> > elemTags, elemNodeTags;
                    gmsh::model::mesh::getElements(elemTypes, elemTags, elemNodeTags, p.first, p.second);
                    if (std::all_of(elemTypes.begin(), elemTypes.end(), [](int &e_type){return e_type == TYP_ELEM_TRIANGLE;} ))
                        { // all elements are triangles
                        int idx = mySets.findFacetteRegionIdx(n);
                        const int nbNod = node.size();
                        for(unsigned int i=0;i<elemNodeTags[0].size();i+=3)
                            {
                            int i0 = elemNodeTags[0][i];
                            int i1 = elemNodeTags[0][i+1];
                            int i2 = elemNodeTags[0][i+2];
                            if (idx > -1) fac.push_back(Facette::Fac(node, nbNod, idx, {i0,i1,i2}));
                            }
                        }
                    else
                        { std::cout << "Fatal error: a surface mesh entity contains other elements than triangles.\n"; exit(1); }
                    }
                }
            }
        } );// end for_each on entities
    }

void mesh::readMesh(Settings const &mySets)
    {
    if (mySets.verbose)
        {
        std::cout << "Reading mesh file " << mySets.getPbName() << ":\n";
        }
    std::string symb;
    std::ifstream msh(mySets.getPbName());

    on_fail_msg_error(msh, "cannot open file " + mySets.getPbName() );

    msh >> symb;
    if (symb == tags::msh::format)
        {
        std::string mshFormat = "";
        msh >> mshFormat;
        if (mshFormat == tags::msh::version)
            {
            int _tags, reg, TYP,nbRegNames;

            if (tags::lookFor(false, msh, tags::msh::begin_physical_names))
                {
                msh >> nbRegNames;

                while ((msh >> symb) && (symb != tags::msh::end_physical_names) && (!msh.fail()))
                    {
                    std::string name;
                    msh >> _tags >> name;

                    switch (stoi(symb))
                        {
                        case tags::msh::DIM_OBJ_2D:
                            {
                            surfRegNames[_tags] = name.substr(1, name.length() - 2);
                            break;
                            }
                        case tags::msh::DIM_OBJ_3D:
                            {
                            volRegNames[_tags] = name.substr(1, name.length() - 2);
                            break;
                            }
                        default:
                            std::cerr << "unknown type in mesh " << tags::msh::begin_physical_names << std::endl;
                            break;
                        }
                    }
                if (mySets.verbose)
                    {
                    std::cout << "  found " << nbRegNames << " regions:\n";
                    std::map<int, std::string>::iterator it;
                    for (it = surfRegNames.begin(); it != surfRegNames.end(); ++it)
                        {
                        std::cout << "    " << it->first << ": " << it->second << '\n';
                        }
                    for (it = volRegNames.begin(); it != volRegNames.end(); ++it)
                        {
                        std::cout << "    " << it->first << ": " << it->second << '\n';
                        }
                    }
                }
            else
                {
                std::cerr << tags::msh::begin_physical_names <<" undefined." << std::endl;
                msh.clear();
                msh.seekg(std::ios::beg);
                }

            double scale = mySets.getScale();
            int nbNod;

            if (tags::lookFor(true, msh, tags::msh::begin_nodes))
                {
                msh >> nbNod;
                if (nbNod > 0) node.resize(nbNod);

                for (int i = 0; i < nbNod; i++)
                    {
                    double x,y,z;
                    msh >> symb >> x >> y >> z;
                    x *= scale;
                    y *= scale;
                    z *= scale;
                    node[i].p = Eigen::Vector3d(x,y,z);
                    }
                }

            on_fail_msg_error(msh, "error while reading nodes");

            for (auto it = surfRegNames.begin(); it != surfRegNames.end(); ++it)
                {
                s.push_back(Mesh::Surf(node, it->second));
                }

            int nbElem;

            if (tags::lookFor(true, msh, tags::msh::begin_elements))
                {
                msh >> nbElem;
                if (mySets.verbose)
                    {
                    std::cout << "  element count: " << nbElem << '\n';
                    }
                while ((msh >> symb) && (symb != tags::msh::end_elements) && (!msh.fail()))
                    {
                    msh >> TYP >> _tags >> reg;
                    for (int i = 1; i < _tags; i++)
                        msh >> symb;
                    switch (TYP)
                        {
                        case tags::msh::TYP_ELEM_TRIANGLE:
                            {
                            int i0, i1, i2;
                            msh >> i0 >> i1 >> i2;
                            if (auto search = surfRegNames.find(reg); search != surfRegNames.end())
                                {  // found named surface
                                for (auto it = s.begin(); it != s.end(); ++it)
                                    {
                                    if (it->getName() == search->second)
                                        {
                                        it->push_back(Mesh::Triangle(node, i0, i1, i2));
                                        }
                                    }

                                int idx = mySets.findFacetteRegionIdx(search->second);
                                if (idx > -1)
                                    fac.push_back(
                                            Facette::Fac(node, nbNod, idx, {i0, i1,
                                                         i2}));  // we only want to store in fac
                                                                // vector the facettes for micromag
                                                                // problem, nothing for bc for stt
                                }
                            else
                                {
                                std::cout << "mesh reading error : unnamed surface region"
                                          << std::endl;
                                }  // unnamed surface

                            break;
                            }
                        case tags::msh::TYP_ELEM_TETRAEDRON:
                            {
                            int i0,i1,i2,i3;
                            msh >> i0 >> i1 >> i2 >> i3;
                            if (auto search = volRegNames.find(reg); search != volRegNames.end())
                                {  // found named volume
                                int idx = mySets.findTetraRegionIdx(search->second);
                                if (idx > -1) tet.push_back(Tetra::Tet(node, idx, {i0,i1,i2,i3}));
                                }
                            else
                                {
                                std::cout << "mesh reading error : unnamed volume region"
                                          << std::endl;
                                }
                            break;
                            }
                        default: std::cerr << "unknown object type in mesh\n"; break;
                        }
                    }
                }
            on_fail_msg_error(msh, "error while reading elements");
            }
        else
            {
            std::cout << "mesh file format " << mshFormat << " not supported." << std::endl;
            SYSTEM_ERROR;
            }
        }

    if (!msh.fail())
        {
        if (mySets.verbose) std::cout << "  closing file\n";
        msh.close();
        }
    else
        {
        std::cout << "error before closing mesh." << std::endl;
        SYSTEM_ERROR;
        }

    for (unsigned int i = 0; i < tet.size(); i++)
        {
        tet[i].idx = i;
        }
    }

double mesh::readSol(bool VERBOSE, const std::string fileName)
    {
    double t(0);
    std::ifstream fin(fileName, std::ifstream::in);

    on_fail_msg_error(fin, "cannot open file " + fileName);

    std::string str;
    bool flag_t = false;
    unsigned long idx;

    while (fin.peek() == '#' || fin.peek() == '\n')
        {
        getline(fin, str);
        idx = str.find(tags::sol::time);
        if (idx != std::string::npos)
            {  // found tag "## time:"
            t = stod(str.substr(idx + tags::sol::time.length() ));
            flag_t = true;
            }
        }

    if (!flag_t)
        {
        std::cerr << "error: no ## time: tag in input .sol file " << fileName << std::endl;
        SYSTEM_ERROR;
        }

    if (VERBOSE)
        {
        std::cout << ".sol file: " << fileName << " @ time t = " << t << std::endl;
        }

    for (unsigned int i = 0; i < node.size(); i++)
        {
        double mx,my,mz;
        unsigned int i_;
        int node_idx = node_index[i];
        fin >> i_ >> mx >> my >> mz >> node[node_idx].d[Nodes::NEXT].phi;
        node[node_idx].d[Nodes::NEXT].u = Eigen::Vector3d(mx,my,mz);
        if (i != i_)
            {
            std::cerr << "error: mesh node index mismatch between mesh and input .sol file"
                      << std::endl;
            fin.close();
            SYSTEM_ERROR;
            }
        }
    fin.close();

    return t;
    }

    }  // namespace Mesh
