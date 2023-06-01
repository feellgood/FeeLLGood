#!/usr/bin/env python3

import numpy as np
import gmsh

def coucou():
    print("coucou toi!")

class mesh(object):
    def __init__ (self,fileName):
        
        self.Tet = []
        self.Names = []
        
        gmsh.initialize()
        gmsh.option.setNumber("General.Terminal",False) # to silent gmsh
        gmsh.open(fileName)
        
        self.readNodes()
        self.readTetra()
        
        gmsh.clear()
        gmsh.finalize()
        
    def readNodes(self):
        if gmsh.model.getDimension() != 3:
            print('Error : convert2vtk does not support Model ' + gmsh.model.getCurrent() + ' : it is not 3D')
        
        # Get all the mesh nodes:
        nodeTags, nodeCoords, nodeParams = gmsh.model.mesh.getNodes()
        self.Nodes = np.empty([len(nodeTags),3],dtype=float)

        for n in range(0,len(nodeTags)):
            idx = int(nodeTags[n])-1
            self.Nodes[idx][0] = nodeCoords[3*n]
            self.Nodes[idx][1] = nodeCoords[3*n + 1]
            self.Nodes[idx][2] = nodeCoords[3*n + 2]
    
    def readTetra(self):    
        entities_3D = gmsh.model.getEntities(dim = 3)
        for e in entities_3D:
            dim = 3
            tag = e[1]
            
            # Get the mesh elements for the entity (dim, tag):
            elemTypes, elemTags, elemNodeTags = gmsh.model.mesh.getElements(dim, tag)
            
            if len(elemTypes) != 1:
                print("Error : convert2vtk only handles mesh with first order tetrahedron elements")
                exit(1)
    
            for n in range(0,len(elemTags[0])):
                idx = int(elemTags[0][n])-1 
                i0 = int(elemNodeTags[0][4*n])
                i1 = int(elemNodeTags[0][4*n+1])
                i2 = int(elemNodeTags[0][4*n+2])
                i3 = int(elemNodeTags[0][4*n+3])
                self.Tet.append([ i0-1 , i1-1, i2-1, i3-1 ])
            # * Does the entity belong to physical groups?
            physicalTags = gmsh.model.getPhysicalGroupsForEntity(dim, tag)
            if len(physicalTags):
                s = ''
                for p in physicalTags:
                    self.Names.append( gmsh.model.getPhysicalName(dim, p))

    def infos(self):
        print("nb Nodes: ",len(self.Nodes),"\tnb Tetrahedrons: ",len(self.Tet), "\tphysical names: ", self.Names)

