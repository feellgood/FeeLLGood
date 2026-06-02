#define BOOST_TEST_MODULE readMeshTest

#include <boost/test/unit_test.hpp>

#include <random>
#include <iostream>
#include <gmsh.h>

#include "tags.h"
#include "mesh.h"
#include "ut_tools.h"
#include "ut_config.h"

// see https://gitlab.onelab.info/gmsh/gmsh/blob/gmsh_4_10_5/tutorials/c++/x1.cpp for more about reading .msh files

BOOST_AUTO_TEST_SUITE(ut_readMesh)

BOOST_AUTO_TEST_CASE(read_ellipsoid)
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
    BOOST_CHECK(nodeT.size() == 167);

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
            BOOST_CHECK(elemTags.size() == 499);
            BOOST_CHECK(nodeTags.size() == 1996);
            }
        else if (pGroup.first == DIM_OBJ_2D)
            {
            std::vector<std::size_t> elemTags, nodeTags;
            gmsh::model::mesh::getElementsByType(TYP_ELEM_TRIANGLE,elemTags, nodeTags );
            std::cout << "2D:\n\telemTags.size= " << elemTags.size() << std::endl;
            std::cout << "\tnodeTags.size= " << nodeTags.size() << std::endl;
            BOOST_CHECK(elemTags.size() == 274);
            BOOST_CHECK(nodeTags.size() == 822);
            }
        } );
    gmsh::clear();
    gmsh::finalize();
    }

BOOST_AUTO_TEST_CASE(checkTriangles)
    {
    // Init basic settings
    Settings shortSet;
    shortSet.paramTetra.emplace_back().regName = "whole_volume";
    shortSet.paramTriangle.emplace_back().regName = "whole_surface";

    std::cout << "Checking the cuboid with a duplicated surface region triangle...\n";
    shortSet.setPbName("meshes/bad_cuboid_1.msh");
    Mesh::mesh badCuboid1(shortSet, true);
    std::stringstream logOutput;
    std::streambuf *buf = std::cerr.rdbuf(logOutput.rdbuf());   // Capture the stderr log output
    bool isCorrect = badCuboid1.checkTriangles();
    BOOST_CHECK(!isCorrect);
    BOOST_CHECK(logOutput.str() == "Error: bad mesh. 2 instances of the same surface triangle have been found\n");
    std::cout << "Done checking the cuboid with a duplicated surface region triangle\n";

    std::cout << "Checking the cuboid with a triangle which isn't in a tetrahedron...\n";
    shortSet.setPbName("meshes/bad_cuboid_2.msh");
    Mesh::mesh badCuboid2(shortSet, true);
    logOutput.str("");
    buf = std::cerr.rdbuf(logOutput.rdbuf());   // Capture the stderr log output
    isCorrect = badCuboid2.checkTriangles();
    BOOST_CHECK(!isCorrect);
    BOOST_CHECK(logOutput.str() == "Error: bad mesh. A triangle which belongs "
                                   "to 1 surface region and no tetrahedron has been found\n");
    std::cout << "Done checking the cuboid with a triangle which isn't in a tetrahedron\n";

    std::cout << "Checking the cuboid with 3 identical triangles in the same tetrahedron...\n";
    shortSet.setPbName("meshes/bad_cuboid_3.msh");
    Mesh::mesh badCuboid3(shortSet, true);
    logOutput.str("");
    buf = std::cerr.rdbuf(logOutput.rdbuf());   // Capture the stderr log output
    isCorrect = badCuboid3.checkTriangles();
    BOOST_CHECK(!isCorrect);
    BOOST_CHECK(logOutput.str() == "Error: bad mesh. A triangular face shared by 3 tetrahedrons has been found\n");
    std::cout << "Done checking the cuboid with 3 identical triangles in the same tetrahedron\n";

    std::cout << "Checking the cuboid with an internal surface region triangle...\n";
    shortSet.setPbName("meshes/bad_cuboid_4.msh");
    Mesh::mesh badCuboid4(shortSet, true);
    logOutput.str("");
    buf = std::cerr.rdbuf(logOutput.rdbuf());   // Capture the stderr log output
    isCorrect = badCuboid4.checkTriangles();
    BOOST_CHECK(!isCorrect);
    BOOST_CHECK(logOutput.str() == "Error: bad mesh. An internal triangle has been "
                                   "found in the surface region 1\n");
    std::cout << "Done checking the cuboid with an internal surface region triangle\n";

    std::cout << "Done testing checkTriangles()\n";
    }

BOOST_AUTO_TEST_SUITE_END()
