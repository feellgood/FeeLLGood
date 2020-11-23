# python -m pip install -U pip
# sudo python -m pip install -U matplotlib

from math import sqrt
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import axes3d
from mpl_toolkits.mplot3d.art3d import Poly3DCollection

import numpy as np

nbX = 3
nbY = 3
nbZ = 1

x_min = -1
x_max = 1
dx = (x_max-x_min)/(nbX)
y_min = -2
y_max = 2
dy = (y_max-y_min)/(nbY)
z_min = -0.6
z_max = 0.6
dz = z_max-z_min

class Cuboid(object):
    def __init__ (self):
        idx = 1
        self.pts = []
        for i in range(0,nbX+1):
            for j in range(0,nbY+1):
                for k in range(0,nbZ+1):
                    self.pts.append([idx,x_min+(i*dx),y_min+(j*dy),z_min+(k*dz)])
                    idx += 1
        self.Tet = []
        self.Fac = []
        for i in range(0,nbX):
            for j in range(0,nbY):
                for k in range(0,nbZ):
                    A = self.calc_idx(i,j,k)
                    B = self.calc_idx(i+1,j,k)
                    C = self.calc_idx(i+1,j+1,k)
                    D = self.calc_idx(i,j+1,k)
                    E = self.calc_idx(i,j,k+1)
                    F = self.calc_idx(i+1,j,k+1)
                    G = self.calc_idx(i+1,j+1,k+1)
                    H = self.calc_idx(i,j+1,k+1)
                    self.Tet.append([A,B,D,H])
                    self.Tet.append([A,B,E,H])
                    self.Tet.append([E,F,H,B])
                    self.Tet.append([B,C,D,H])
                    self.Tet.append([B,C,G,H])
                    self.Tet.append([F,B,G,H])

                    if (nbZ == 1): 
                        self.Fac.append([A,B,D])
                        self.Fac.append([B,C,D])
                        self.Fac.append([E,F,H])
                        self.Fac.append([F,G,H])
                    elif (k == 0):
                        self.Fac.append([A,B,D])
                        self.Fac.append([B,C,D])
                    elif (k == (nbZ-1)):
                        self.Fac.append([E,F,H])
                        self.Fac.append([F,G,H])

                    if (nbX == 1): 
                        self.Fac.append([A,D,H])
                        self.Fac.append([A,E,H])
                        self.Fac.append([B,C,G])
                        self.Fac.append([B,F,H])
                    elif (i == 0):
                        self.Fac.append([A,D,H])
                        self.Fac.append([A,E,H])
                    elif (i == (nbX-1)):
                        self.Fac.append([B,C,G])
                        self.Fac.append([B,F,H])

                    if (nbY == 1): 
                        self.Fac.append([A,B,E])
                        self.Fac.append([B,E,F])
                        self.Fac.append([C,G,H])
                        self.Fac.append([H,D,C])
                    elif (j == 0):
                        self.Fac.append([A,B,E])
                        self.Fac.append([B,E,F])
                    elif (j == (nbY-1)):
                        self.Fac.append([C,G,H])
                        self.Fac.append([H,D,C])

    def stringTet(self,idx):
        return str(1+self.Tet[idx][0])+"\t"+str(1+self.Tet[idx][1])+"\t"+str(1+self.Tet[idx][2]) +"\t"+str(1+self.Tet[idx][3])

    def stringFac(self,idx):
        return str(1+self.Fac[idx][0])+"\t"+str(1+self.Fac[idx][1])+"\t"+str(1+self.Fac[idx][2])

    def calc_idx(self,i,j,k):
        return( (nbZ+1)*(nbY+1)*i +  (nbZ+1)*j + k)

    def make(self,meshFileName,volRegionName,surfRegionName):
        meshFile = open(meshFileName,'w')
        meshFile.write("$MeshFormat\n2.2\t0\t8\n$EndMeshFormat\n$Nodes\n")
        meshFile.write(str(len(self.pts))+"\n")
        for i in range(0,len(self.pts)):
            meshFile.write( str(self.pts[i][0])+"\t"+str(self.pts[i][1])+"\t"+str(self.pts[i][2])+"\t"+str(self.pts[i][3]) + "\n" )
        meshFile.write("$EndNodes\n$Elements\n")
        meshFile.write(str(len(self.Tet)+len(self.Fac))+"\n")
        idx=1
        for i in range(0,len(self.Tet)):
            s = self.stringTet(i)
            meshFile.write( str(idx)+"\t4\t2\t" + volRegionName + "\t1\t" + s + "\n" )
            idx += 1
        for i in range(0,len(self.Fac)):
            s = self.stringFac(i)
            meshFile.write( str(idx)+"\t2\t2\t" + surfRegionName + "\t1\t"+ s + "\n" )
            idx += 1
        meshFile.write("$EndElements\n")
        meshFile.close()
        print("mesh file " + meshFileName + " generated.")

rectangle = Cuboid()
rectangle.make("rectangle.msh","300","200")


