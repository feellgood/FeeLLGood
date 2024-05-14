#define BOOST_TEST_MODULE readMeshTest

#include <boost/test/unit_test.hpp>

#include <random>
#include <iostream>
#include <gmsh.h>

#include "tags.h"
#include "ut_tools.h"
#include "ut_config.h"

// see https://gitlab.onelab.info/gmsh/gmsh/blob/gmsh_4_10_5/tutorials/c++/x1.cpp for more about reading .msh files

BOOST_AUTO_TEST_SUITE(ut_readMesh)

BOOST_AUTO_TEST_CASE(readMesh)
    {
    using namespace tags::msh;
    gmsh::initialize();

    std::string fileName("../examples/ellipsoid.msh");
    gmsh::open(fileName);

    std::string name;
    gmsh::model::getCurrent(name);
    int dim = gmsh::model::getDimension();
    std::cout << "Model " << name << " (" << dim << "D)\n";
    BOOST_CHECK(name == "ellipsoid");
    BOOST_CHECK(dim == DIM_OBJ_3D);

    std::vector<std::size_t> nodeT;
    std::vector<double> nodeC, nodeP;
    gmsh::model::mesh::getNodes(nodeT, nodeC, nodeP, -1,-1); // all nodes of the mesh
    BOOST_CHECK(nodeT.size() == 163);

    std::vector<int> elemTypes;
    gmsh::model::mesh::getElementTypes(elemTypes, DIM_OBJ_3D, -1);
    BOOST_CHECK(elemTypes.size() == 1); // only one type of element allowed
    BOOST_CHECK(elemTypes[0] == TYP_ELEM_TETRAEDRON); // only Tetrahedrons

    gmsh::model::mesh::getElementTypes(elemTypes, DIM_OBJ_2D, -1);
    BOOST_CHECK(elemTypes.size() == 1); // only one type of element allowed
    BOOST_CHECK(elemTypes[0] == TYP_ELEM_TRIANGLE); // only Triangles

    std::vector<std::pair<int, int> > physGroups;
    gmsh::model::getPhysicalGroups(physGroups,-1); // -1 for all
    std::for_each(physGroups.begin(),physGroups.end(),[]( std::pair<int, int> &pGroup)
        {
        std::string physGroupName;
        gmsh::model::getPhysicalName(pGroup.first, pGroup.second, physGroupName);
        std::cout << pGroup.first << " ; "<< pGroup.second <<  " Name= " << physGroupName << std::endl;
        if (pGroup.first == DIM_OBJ_3D)
            {
            std::vector<std::size_t> elemTags, nodeTags;
            gmsh::model::mesh::getElementsByType(TYP_ELEM_TETRAEDRON,elemTags, nodeTags );
            std::cout << "3D:\n\telemTags.size= " << elemTags.size() << std::endl;
            std::cout << "\tnodeTags.size= " << nodeTags.size() << std::endl;
            BOOST_CHECK(elemTags.size() == 320);
            BOOST_CHECK(nodeTags.size() == 1280);
            }
        else if (pGroup.first == DIM_OBJ_2D)
            {
            std::vector<std::size_t> elemTags, nodeTags;
            gmsh::model::mesh::getElementsByType(TYP_ELEM_TRIANGLE,elemTags, nodeTags );
            std::cout << "2D:\n\telemTags.size= " << elemTags.size() << std::endl;
            std::cout << "\tnodeTags.size= " << nodeTags.size() << std::endl;
            BOOST_CHECK(elemTags.size() == 320);
            BOOST_CHECK(nodeTags.size() == 960);
            }
        } );
    gmsh::clear();
    gmsh::finalize();
    }

BOOST_AUTO_TEST_SUITE_END()

