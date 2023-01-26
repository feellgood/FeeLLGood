#!/usr/bin/env python3

import math
import gmsh

import numpy as np

class Cylinder(object):
    def __init__ (self,radius,thickness,mesh_size,surfName,volName):
        """ 
            geometrical cylinder is zero centered, with radius r and length t along (Oz), build by extrusion
            gmsh file format 2.2 is used to write the mesh text file  
        """
        
        self.r = radius
        self.t = thickness
        self.msh_s = mesh_size
        self.surfName = surfName
        self.volName = volName
        self.withExtraSurf = False

    def addEdgeSurf(self,name1,name2):
        """ optional surfaces : name1 will refer to the base surface, name2 will refer to the translated name1 surface from extrusion """
        self.n1 = name1
        self.n2 = name2
        self.withExtraSurf = True
    
    def make(self,meshFileName):
        """ write cylinder mesh file in gmsh 2.2 text format """
        
        gmsh.initialize()
        gmsh.option.setNumber("General.Terminal",False) # to silent gmsh
        gmsh.model.add("cyl")
        
        p_origin = gmsh.model.geo.addPoint(0,0,-0.5*self.t,self.msh_s)

        start_pt = gmsh.model.geo.addPoint(self.r,0,-0.5*self.t,self.msh_s)
        interim_pt = gmsh.model.geo.addPoint(0,self.r,-0.5*self.t,self.msh_s)
        interim1_pt = gmsh.model.geo.addPoint(-self.r,0,-0.5*self.t,self.msh_s)
        interim2_pt = gmsh.model.geo.addPoint(0,-self.r,-0.5*self.t,self.msh_s)
        end_pt = gmsh.model.geo.addPoint(self.r,0,-0.5*self.t,self.msh_s)

        big_circle = gmsh.model.geo.addCircleArc(start_pt,p_origin,interim_pt)
        big_circle1 = gmsh.model.geo.addCircleArc(interim_pt,p_origin,interim1_pt)
        big_circle2 = gmsh.model.geo.addCircleArc(interim1_pt,p_origin,interim2_pt)
        big_circle3 = gmsh.model.geo.addCircleArc(interim2_pt,p_origin,end_pt)

        curvedLoop = gmsh.model.geo.addCurveLoop([big_circle,big_circle1,big_circle2,big_circle3]) # curvedLoop is an index (integer)

        surf = gmsh.model.geo.addPlaneSurface([curvedLoop]) # surf is an index (integer)

        out = gmsh.model.geo.extrude([(2,surf)],0,0,self.t) # 2 is the dimension of the object refered by index surf

        gmsh.model.geo.synchronize() # we have to sync before calling addPhysicalGroup

        surface_tag = 200
        gmsh.model.addPhysicalGroup(2,[surf, out[0][1], out[2][1], out[3][1], out[4][1], out[5][1]],surface_tag)
        gmsh.model.setPhysicalName(2,surface_tag,self.surfName)

        if self.withExtraSurf :
                surface_tag_left = 201
                gmsh.model.addPhysicalGroup(2,[surf],surface_tag_left)
                gmsh.model.setPhysicalName(2,surface_tag_left,self.n1)

                surface_tag_right = 202
                gmsh.model.addPhysicalGroup(2,[ out[0][1] ],surface_tag_right)
                gmsh.model.setPhysicalName(2,surface_tag_right,self.n2)

        volume_tag = 300
        gmsh.model.addPhysicalGroup(3,[out[1][1]],volume_tag)
        gmsh.model.setPhysicalName(3,volume_tag,self.volName)

        gmsh.model.geo.synchronize() # we have to synchronize before the call to 'generate' to build the mesh
        gmsh.model.mesh.generate(3) # 3 is the dimension of the mesh

        gmsh.option.set_number("Mesh.MshFileVersion", 2.2) # to force mesh file format to 2.2

        gmsh.write(meshFileName)

        nbNodes = gmsh.option.get_number("Mesh.NbNodes")
        nbTriangles = gmsh.option.get_number("Mesh.NbTriangles")
        nbTetra = gmsh.option.get_number("Mesh.NbTetrahedra")

# uncomment next line to see a graphic rendering of the mesh
        #gmsh.fltk.run()
                
        gmsh.finalize()
        
        print("mesh file " + meshFileName + " generated, nb nodes =" + str(nbNodes) + " , nb tetra =" +str(nbTetra) )

