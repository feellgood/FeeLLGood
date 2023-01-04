import numpy as np

class Cuboid(object):
    def __init__ (self,pt_min,pt_max,nbX,nbY,nbZ):
        """ 
            constructor, inits nodes lists, and builds tetrahedrons and outer surface mesh of the cuboid in the volume defined by pt_min,pt_max
            Tet is a table containing 4 indices refering to the list of points
            Fac is a table containing 3 indices refering to the list of points
            index shift to obey gmsh file format 2.2 is performed while writing the file, except first column indices of points table  
        """
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
                    self.pts.append([idx,pt_min[0]+(i*dx),pt_min[1]+(j*dy),pt_min[2]+(k*dz)])
                    idx += 1
        self.Tet = []
        self.Fac = []
        self.subSurf = []
        self.subSurfRegName = []
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

    def transformNodes(self,func):
        """
        apply func on pts array, func input and output are 3D vectors
        """
        for i in range(0,len(self.pts)):
            # pts[i][0] is an index, so we do not use map function here
            x = self.pts[i][1]
            y = self.pts[i][2]
            z = self.pts[i][3]
            transformed_pt = []
            transformed_pt = func([x,y,z])
            self.pts[i][1] = transformed_pt[0]
            self.pts[i][2] = transformed_pt[1]
            self.pts[i][3] = transformed_pt[2]
            
    def stringTet(self,idx):
        return str(1+self.Tet[idx][0])+"\t"+str(1+self.Tet[idx][1])+"\t"+str(1+self.Tet[idx][2]) +"\t"+str(1+self.Tet[idx][3])

    def stringFac(self,table,idx):
        return str(1+table[idx][0])+"\t"+str(1+table[idx][1])+"\t"+str(1+table[idx][2])
    
    def calc_idx(self,i,j,k):
        return( (self.nbZ+1)*(self.nbY+1)*i +  (self.nbZ+1)*j + k)

    def calc_nb_elements(self):
        nbElem = len(self.Tet) + len(self.Fac)
        for i in range(0,len(self.subSurf)):
            nbElem += len(self.subSurf[i])
        return(nbElem)
        
    def add_sub_surface(self,subSurfRegName,func):
        tempo_subSurf = []
        for i in range(0,len(self.Fac)):
            f = self.Fac[i]
            idxA = f[0]
            idxB = f[1]
            idxC = f[2]
            A = self.pts[idxA]
            B = self.pts[idxB]
            C = self.pts[idxC]
            
            if(func(A[1],A[2],A[3]) and func(B[1],B[2],B[3]) and func(C[1],C[2],C[3]) ):
                tempo_subSurf.append( [idxA,idxB,idxC] )
        self.subSurf.append( tempo_subSurf )
        self.subSurfRegName.append(subSurfRegName)

    
    def make(self,meshFileName,volRegionName,surfRegionName):
        """ write mesh file in gmsh 2.2 format """
        meshFile = open(meshFileName,'w')
        meshFile.write("$MeshFormat\n2.2\t0\t8\n$EndMeshFormat\n")
        nbTotNames = 2 + len(self.subSurfRegName)
        meshFile.write("$PhysicalNames\n" + str(nbTotNames) + "\n")
        surfRegionTag = 200
        if (len(self.subSurfRegName) < 99):
            volRegionTag = 300 
        else:
            volRegionTag = 300 + len(self.subSurfRegName)
        meshFile.write("2 " + str(surfRegionTag) + ' \"' + surfRegionName + '\"\n')
        meshFile.write("3 " + str(volRegionTag) + ' \"' + volRegionName + '\"\n')
        
        for i in range(0,len(self.subSurfRegName)):
            meshFile.write( "2 " + str(surfRegionTag + 1 + i) + ' \"' + self.subSurfRegName[i] +'\"\n'  )
        meshFile.write("$EndPhysicalNames\n")
        meshFile.write("$Nodes\n" + str(len(self.pts))+"\n")
        for i in range(0,len(self.pts)):
            meshFile.write( str(self.pts[i][0])+"\t"+str(self.pts[i][1])+"\t"+str(self.pts[i][2])+"\t"+str(self.pts[i][3]) + "\n" )
        meshFile.write("$EndNodes\n$Elements\n")
        nbElem = self.calc_nb_elements()
        meshFile.write(str(nbElem)+"\n")
        idx=1
        for i in range(0,len(self.Tet)):
            s = self.stringTet(i)
            meshFile.write( str(idx)+"\t4\t2\t" + str(volRegionTag) + "\t1\t" + s + "\n" )
            idx += 1
        for i in range(0,len(self.Fac)):
            s = self.stringFac(self.Fac,i)
            meshFile.write( str(idx)+"\t2\t2\t" + str(surfRegionTag) + "\t1\t"+ s + "\n" )
            idx += 1
        for i in range(0,len(self.subSurf)):
            for j in range(0,len(self.subSurf[i])):
                s = self.stringFac(self.subSurf[i],j)
                meshFile.write( str(idx)+"\t2\t2\t" + self.subSurfRegName[i] + "\t1\t"+ s + "\n" )
                idx += 1
        meshFile.write("$EndElements\n")
        meshFile.close()
        print("mesh file " + meshFileName + " generated, nb nodes =" + str(len(self.pts)) + " , nb tetra =" +str(len(self.Tet)) )

