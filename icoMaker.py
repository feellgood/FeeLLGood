# python -m pip install -U pip
# sudo python -m pip install -U matplotlib

from math import sqrt
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import axes3d
from mpl_toolkits.mplot3d.art3d import Poly3DCollection

import numpy as np

# this script contains two functions
# makeIco builds a mesh of an icosahedron generated from scrap, each tetrahedron is a face of the icosahedron linked to its center.
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

# the following method is only to have a 3D view of the icosahedron, useless if you do not want to see it
    def show(self):
        fig = plt.figure()
        ax = fig.add_subplot(111,projection='3d')

        ax.set_xlabel('X')
        ax.set_xlim3d(-2, 2)
        ax.set_ylabel('Y')
        ax.set_ylim3d(-2, 2)
        ax.set_zlabel('Z')
        ax.set_zlim3d(-2, 2)

        r1 = [ self.pts[0], self.pts[2], self.pts[3], self.pts[1] ] 
        r2 = [ self.pts[4], self.pts[5], self.pts[7], self.pts[6] ]
        r3 = [ self.pts[8], self.pts[10], self.pts[11], self.pts[9] ]

        triangle3D = []
        for i in range(0,len(self.f_idx)):
            triangle3D.append( [self.pts[ self.f_idx[i][0] ],self.pts[ self.f_idx[i][1] ],self.pts[ self.f_idx[i][2] ]] )

        triangle3D.append(r1)
        triangle3D.append(r2)
        triangle3D.append(r3)

        fc = []
        for i in range(1,len(self.f_idx)+1):
            fc.append( (i/20,0,0,0.2) )

        fc.append( (1,1,0,0.4) )
        fc.append( (0,1,0,0.4) )
        fc.append( (0,0,1,0.4) )

        ax.add_collection3d( Poly3DCollection( triangle3D, facecolors=fc, linewidths=1))
        plt.show()
