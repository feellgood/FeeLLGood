#!/usr/bin/env python3
# python -m pip install -U pip

from math import sqrt

import numpy as np

# Ico class
# make builds a mesh of an icosahedron format2.2, each tetrahedron is built with a outer face of the icosahedron linked to the barycenter of the icosahedron (zero). The mesh file is directly written in the output textfile, without calling gmsh module.
# https://fr.wikipedia.org/wiki/Icosa%C3%A8dre

class Ico(object):
    def __init__ (self):
        self.phi = (1+sqrt(5))/2
        self.pts = np.array([[self.phi,1,0], [self.phi,-1,0], [-self.phi,1,0], [-self.phi,-1,0],[1,0,self.phi], [1,0,-self.phi], [-1,0,self.phi], [-1,0,-self.phi],[0,self.phi,1], [0,self.phi,-1], [0,-self.phi,1], [0,-self.phi,-1], [0,0,0]] )

        self.f_idx = [ [0,8,4], [0,4,1], [0,1,5], [0,5,9], [0,9,8], [8,9,2], [8,2,6], [8,6,4], [6,10,4], [6,3,10], [6,2,3], [10,3,11], [10,11,1], [10,1,4], [1,11,5], [11,3,7], [11,7,5], [3,2,7], [7,9,2],  [5,7,9] ]

    def make(self,meshFileName,volRegionName,surfRegionName):
        meshFile = open(meshFileName,'w')
        meshFile.write("$MeshFormat\n2.2\t0\t8\n$EndMeshFormat\n$Nodes\n")
        meshFile.write(str(len(self.pts))+"\n")
        for i in range(0,len(self.pts)):
            meshFile.write( str(i+1)+"\t"+str(self.pts[i][0])+"\t"+str(self.pts[i][1])+"\t"+str(self.pts[i][2]) + "\n" )
        meshFile.write("$EndNodes\n$Elements\n")
        meshFile.write(str(2*len(self.f_idx))+"\n")
        for i in range(0,len(self.f_idx)):
            meshFile.write( str(i+1)+"\t4\t2\t" + volRegionName + "\t1\t13\t"+str(1+self.f_idx[i][0])+"\t"+str(1+self.f_idx[i][1])+"\t"+str(1+self.f_idx[i][2]) + "\n" )
        for i in range(0,len(self.f_idx)):
            meshFile.write( str(len(self.f_idx)+i+1)+"\t2\t2\t" + surfRegionName + "\t1\t"+str(1+self.f_idx[i][0])+"\t"+str(1+self.f_idx[i][1])+"\t"+str(1+self.f_idx[i][2]) + "\n" )

        meshFile.write("$EndElements\n")
        meshFile.close()
        print("mesh file " + meshFileName + " generated.")


