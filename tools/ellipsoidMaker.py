#!/usr/bin/env python3

import math
import gmsh

import numpy as np

class Ellipsoid(object):
    def __init__ (self,r1,r2,mesh_size,surfName,volName):
        """ 
            geometrical ellipsoid is zero centered, with radius r1 in the plane (Oxy) and r2 along (Oz)
            gmsh file format 2.2 is used to write the mesh text file  
        """
        
        self.r1 = r1
        self.r2 = r2
        self.msh_s = mesh_size
        self.surfName = surfName
        self.volName = volName
    
    def make(self,meshFileName):
        """ write ellpipsoid mesh file in gmsh 2.2 text format """
        
        gmsh.initialize()

        gmsh.model.add("ellipsoid")
        sph = gmsh.model.occ.addSphere(0,0,0,self.r1)
        gmsh.model.occ.dilate([(3,sph)],0,0,0,1,1,self.r2)
        gmsh.model.occ.synchronize() # we have to sync before calling addPhysicalGroup

        surfaces = gmsh.model.getEntities(dim=2)
        surface_tag = 200
        gmsh.model.addPhysicalGroup(surfaces[0][0], [surfaces[0][1]], surface_tag) # first parameter should be 2
        gmsh.model.setPhysicalName(2,surface_tag,self.surfName)

        volumes = gmsh.model.getEntities(dim=3)
        volume_tag = 300
        gmsh.model.addPhysicalGroup(volumes[0][0], [volumes[0][1]], volume_tag) # first parameter should be 3
        gmsh.model.setPhysicalName(3,volume_tag,self.volName)

        gmsh.model.occ.synchronize() # we have to synchronize before the call to 'generate' to build the mesh
        
        eps = 1e-3
        xmin = -self.r1 - eps
        ymin = -self.r1 - eps
        zmin = -self.r2 - eps
        xmax = self.r1 + eps
        ymax = self.r1 + eps
        zmax = self.r2 + eps
        
        objects = gmsh.model.getEntitiesInBoundingBox(xmin,ymin,zmin,xmax,ymax,zmax, dim = 0)
        gmsh.model.mesh.setSize(objects,self.msh_s)
        
        gmsh.model.mesh.generate(dim = 3)

        gmsh.option.set_number("Mesh.MshFileVersion", 2.2) # to force mesh file format to 2.2

        gmsh.write(meshFileName)

        nbNodes = gmsh.option.get_number("Mesh.NbNodes")
        nbTriangles = gmsh.option.get_number("Mesh.NbTriangles")
        nbTetra = gmsh.option.get_number("Mesh.NbTetrahedra")

# uncomment next line to see a graphic rendering of the mesh
        #gmsh.fltk.run()
                
        gmsh.finalize()
        
        print("mesh file " + meshFileName + " generated, nb nodes =" + str(nbNodes) + " , nb tetra =" +str(nbTetra) )

