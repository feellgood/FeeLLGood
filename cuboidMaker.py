# python -m pip install -U pip
# sudo python -m pip install -U matplotlib

import numpy as np

class Cuboid(object):
    def __init__ (self,pt_min,pt_max,nbX,nbY,nbZ):
        idx = 1
        self.nbX = nbX
        self.nbY = nbY
        self.nbZ = nbZ
        dx = (pt_max[0]-pt_min[0])/self.nbX
        dy = (pt_max[1]-pt_min[1])/self.nbY
        dz = (pt_max[2]-pt_min[2])/self.nbZ
        self.pts = []
        for i in range(0,self.nbX+1):
            for j in range(0,self.nbY+1):
                for k in range(0,self.nbZ+1):
                    self.pts.append([idx,pt_min[0]+(i*dx),pt_min[1]+(j*dy),pt_min[1]+(k*dz)])
                    idx += 1
        self.Tet = []
        self.Fac = []
        for i in range(0,self.nbX):
            for j in range(0,self.nbY):
                for k in range(0,self.nbZ):
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

                    if (self.nbZ == 1): 
                        self.Fac.append([A,B,D])
                        self.Fac.append([B,C,D])
                        self.Fac.append([E,F,H])
                        self.Fac.append([F,G,H])
                    elif (k == 0):
                        self.Fac.append([A,B,D])
                        self.Fac.append([B,C,D])
                    elif (k == (self.nbZ-1)):
                        self.Fac.append([E,F,H])
                        self.Fac.append([F,G,H])

                    if (self.nbX == 1): 
                        self.Fac.append([A,D,H])
                        self.Fac.append([A,E,H])
                        self.Fac.append([B,C,G])
                        self.Fac.append([B,F,G])
                    elif (i == 0):
                        self.Fac.append([A,D,H])
                        self.Fac.append([A,E,H])
                    elif (i == (self.nbX-1)):
                        self.Fac.append([B,C,G])
                        self.Fac.append([B,F,G])

                    if (self.nbY == 1): 
                        self.Fac.append([A,B,E])
                        self.Fac.append([B,E,F])
                        self.Fac.append([C,G,H])
                        self.Fac.append([H,D,C])
                    elif (j == 0):
                        self.Fac.append([A,B,E])
                        self.Fac.append([B,E,F])
                    elif (j == (self.nbY-1)):
                        self.Fac.append([C,G,H])
                        self.Fac.append([H,D,C])

    def stringTet(self,idx):
        return str(1+self.Tet[idx][0])+"\t"+str(1+self.Tet[idx][1])+"\t"+str(1+self.Tet[idx][2]) +"\t"+str(1+self.Tet[idx][3])

    def stringFac(self,idx):
        return str(1+self.Fac[idx][0])+"\t"+str(1+self.Fac[idx][1])+"\t"+str(1+self.Fac[idx][2])

    def calc_idx(self,i,j,k):
        return( (self.nbZ+1)*(self.nbY+1)*i +  (self.nbZ+1)*j + k)

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

rectangle = Cuboid([-1,-2,-3],[1,2,3],4,4,2)
rectangle.make("rectangle.msh","300","200")


