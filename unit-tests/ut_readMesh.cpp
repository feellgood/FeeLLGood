#define BOOST_TEST_MODULE readMeshTest

#include <boost/test/unit_test.hpp>

#include <random>
#include <iostream>
#include <gmsh.h>

#include "ut_tools.h"
#include "ut_config.h"

// very inspired of https://gitlab.onelab.info/gmsh/gmsh/blob/gmsh_4_10_5/tutorials/c++/x1.cpp

BOOST_AUTO_TEST_SUITE(ut_readMesh)

BOOST_AUTO_TEST_CASE(readMesh)
    {
    gmsh::initialize();
    
    std::string fileName("../examples/ellipsoid.msh");
    gmsh::open(fileName);
    
    std::string name;
    gmsh::model::getCurrent(name);
    int dim = gmsh::model::getDimension();
    std::cout << "Model " << name << " (" << dim << "D)\n";
    BOOST_CHECK(name == "ellipsoid");
    BOOST_CHECK(dim == 3);

    std::vector<std::pair<int, int> > entities;
    gmsh::model::getEntities(entities);

    std::vector<std::size_t> nodeT;
    std::vector<double> nodeC, nodePrms;
    dim = 2;
    int tag=1;
    gmsh::model::mesh::getNodes(nodeT, nodeC, nodePrms, dim, tag);
    BOOST_CHECK(nodeT.size() == 162); // surface has 162 nodes
    
    for(auto e : entities)
        {
        int dim = e.first, tag = e.second;
        std::vector<std::size_t> nodeTags;
        std::vector<double> nodeCoords, nodeParams;
        gmsh::model::mesh::getNodes(nodeTags, nodeCoords, nodeParams, dim, tag);

        std::vector<int> elemTypes;
        std::vector<std::vector<std::size_t> > elemTags, elemNodeTags;
        gmsh::model::mesh::getElements(elemTypes, elemTags, elemNodeTags, dim, tag);

        std::string type;
        gmsh::model::getType(dim, tag, type);
        std::string name;
        gmsh::model::getEntityName(dim, tag, name);
        if(name.size()) name += " ";
        std::cout << "Entity " << name << "(" << dim << "," << tag << ") of type " << type << "\n";

        int numElem = 0;
        for(auto &tags : elemTags) numElem += tags.size();
        std::cout << " - Mesh has " << nodeTags.size() << " nodes and " << numElem << " elements\n";

        std::vector<int> physicalTags;
        gmsh::model::getPhysicalGroupsForEntity(dim, tag, physicalTags);
        if(physicalTags.size())
            {
            std::cout << " - Physical group: ";
            for(auto physTag : physicalTags)
                {
                std::string n;
                gmsh::model::getPhysicalName(dim, physTag, n);
                if(n.size()) n += " ";
                std::cout << n << "(" << dim << ", " << physTag << ") ";
                }
            std::cout << "\n";
            }

        std::vector<int> partitions;
        gmsh::model::getPartitions(dim, tag, partitions);
        if(partitions.size())
            {
            std::cout << " - Partition tags:";
            for(auto part : partitions) std::cout << " " << part;
            int parentDim, parentTag;
            gmsh::model::getParent(dim, tag, parentDim, parentTag);
            std::cout << " - parent entity (" << parentDim << "," << parentTag << ")\n";
            }

        for(auto elemType : elemTypes)
            {
            std::string name;
            int d, order, numv, numpv;
            std::vector<double> param;
            gmsh::model::mesh::getElementProperties(elemType, name, d, order, numv, param, numpv);
            std::cout << " - Element type: " << name << ", order " << order << "\n";
            std::cout << "   with " << numv << " nodes in param coord: (";
            for(auto p : param) std::cout << p << " ";
            std::cout << ")\n";
            }
        }

    gmsh::clear();
    gmsh::finalize();
    }

BOOST_AUTO_TEST_SUITE_END()

