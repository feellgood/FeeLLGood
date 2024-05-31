#include <fstream>
#include <algorithm>
#include <gmsh.h>

#include "feellgoodSettings.h"
#include "mesh.h"
#include "meshUtils.h"
#include "tags.h"

namespace Mesh
    {
void mesh::readMesh(Settings const &mySets)
    {
    gmsh::initialize();
    if (!mySets.verbose)
        { gmsh::option::setNumber("General.Terminal",0); } // to silent gmsh
    checkMeshFile(mySets);
    if (mySets.verbose)
        { std::cout<< "Checking mesh file: done.\n"; }
    readNodes(mySets);
    if (mySets.verbose)
        { std::cout<< "Nodes loaded.\n"; }
    readTetraedrons(mySets);
    for (unsigned int i = 0; i < tet.size(); i++)
        { tet[i].idx = i; }
    if (mySets.verbose)
        { std::cout<< "Tetraedrons built.\n"; }
    readTriangles(mySets);
    if (mySets.verbose)
        { std::cout<< "triangles built.\n"; }
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
    using namespace tags::msh;
    double scale = mySets.getScale();
    std::vector<std::pair<int, int> > physGroups;
    gmsh::model::getPhysicalGroups(physGroups,DIM_OBJ_3D);
    int i(0);
    std::for_each(physGroups.begin(),physGroups.end(),[this,scale,&i](std::pair<int, int> &pGroup)
        {
        std::vector<std::size_t> nodeT;
        std::vector<double> nodeC;
        gmsh::model::mesh::getNodesForPhysicalGroup(DIM_OBJ_3D,pGroup.second,nodeT,nodeC);
        node.resize( node.size() + nodeT.size() );// ouch...
        std::for_each(nodeT.begin(),nodeT.end(),[this,&i,scale,&nodeC](std::size_t &loc_idx)
            {
            long j = loc_idx-1;
            node[i].p = Eigen::Vector3d(nodeC[3*j], nodeC[3*j + 1], nodeC[3*j + 2]) * scale;
            i++;
            });
        });
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
                    if (all_elems_are<TYP_ELEM_TETRAEDRON>(elemTypes))
                        { // all elements are tetrahedrons
                        int idx = mySets.findTetraRegionIdx(n);
                        for(unsigned int i=0;i<elemNodeTags[0].size();i+=SIZE_TETRAEDRON)
                            {
                            int i0 = elemNodeTags[0][i];
                            int i1 = elemNodeTags[0][i+1];
                            int i2 = elemNodeTags[0][i+2];
                            int i3 = elemNodeTags[0][i+3];
                            if (idx > -1)
                                {
                                tet.push_back(Tetra::Tet(node, idx, {i0,i1,i2,i3}));
                                }
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
                    if (all_elems_are<TYP_ELEM_TRIANGLE>(elemTypes))
                        { // all elements are triangles
                        int idx = mySets.findFacetteRegionIdx(n);
                        const int nbNod = node.size();
                        for(unsigned int i=0;i<elemNodeTags[0].size();i+=SIZE_TRIANGLE)
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
