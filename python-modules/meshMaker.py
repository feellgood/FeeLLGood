#!/usr/bin/env python3

from math import sqrt
import numpy as np

def gmsh_init(objName,fileFormat,binary,verbose=True):
    """
    initializations of gmsh used by some classes below. Output file format set to 2.2 or 4.1 and text or binary
    """
    import gmsh
    global gmsh
    gmsh.initialize()
    gmsh.option.setNumber("General.Terminal",verbose) # to silent(or not) gmsh
    file_format_versions = [2.1, 2.2, 3.0, 4.0, 4.1]
    if (fileFormat in file_format_versions):
        gmsh.option.setNumber("Mesh.MshFileVersion", fileFormat)
    gmsh.option.setNumber("Mesh.Binary",binary)
    gmsh.model.add(objName)

def gmsh_mesh_build(name):
    """
    remove duplicate nodes of the mesh (seems to be useless in most situations)
    """
    gmsh.model.setCurrent(name)
    gmsh.model.mesh.removeDuplicateNodes()
    gmsh.model.geo.synchronize() # we have to synchronize before the call to 'generate' to build the mesh
    gmsh.model.mesh.generate(3) # 3 is the dimension of the mesh

def infos(g_opt,name):
    """
    returns a string with total number of nodes of 3D physical entities, number of tetraedrons, and number of triangles corresponding  to written mesh file of model 'name'
    """
    gmsh.model.setCurrent(name)
    e_3D = gmsh.model.getEntities(dim = 3)
    nbNodes = 0
    for o3D in e_3D:
        physicalTags = gmsh.model.getPhysicalGroupsForEntity( 3, o3D[1] )
        for p in physicalTags:
            nodeTags = gmsh.model.mesh.getNodesForPhysicalGroup(3,p)
            nbNodes += int(len(nodeTags[1])/3)
    
    nbTriangles = int(g_opt.getNumber("Mesh.NbTriangles"))
    nbTetra = int(g_opt.getNumber("Mesh.NbTetrahedra"))
    return f"{nbNodes} nodes, {nbTetra} tetrahedra, {nbTriangles} triangles"

class Hexahedron(object):
    def __init__ (self,pt_min,pt_max,mesh_size,surfName,volName,fileFormat=4.1,verbose=False):
        """
        hexahedron defined by pt_min and pt_max, mesh_size average tetrahedron size
        bool optimize is optional optimization of the mesh by relocating of the nodes.
        by default output mesh file is written using text file format 4.1
        """
        self.msh_s = mesh_size
        self.surfName = surfName
        self.volName = volName
        self.p_min = pt_min
        self.p_max = pt_max
        gmsh_init("hex",fileFormat,0,verbose)
        
    def make(self,meshFileName):
        """ write hexahedron mesh file """
        
        p = np.empty(4, dtype=object)
        p[0] = gmsh.model.geo.addPoint(self.p_min[0],self.p_min[1],self.p_min[2],self.msh_s)
        p[1] = gmsh.model.geo.addPoint(self.p_max[0],self.p_min[1],self.p_min[2],self.msh_s)
        p[2] = gmsh.model.geo.addPoint(self.p_max[0],self.p_max[1],self.p_min[2],self.msh_s)
        p[3] = gmsh.model.geo.addPoint(self.p_min[0],self.p_max[1],self.p_min[2],self.msh_s)
        
        l = np.empty(4,dtype=object)
        l[0] = gmsh.model.geo.addLine(p[0],p[1])
        l[1] = gmsh.model.geo.addLine(p[1],p[2])
        l[2] = gmsh.model.geo.addLine(p[2],p[3])
        l[3] = gmsh.model.geo.addLine(p[3],p[0])
        
        closedLoop = gmsh.model.geo.addCurveLoop([l[0],l[1],l[2],l[3]])
        surf = gmsh.model.geo.addPlaneSurface([closedLoop]) # surf is an index (integer)
        t = self.p_max[2] - self.p_min[2] 
        out = gmsh.model.geo.extrude([(2,surf)],0,0,t) # 2 is the dimension of the object refered by index surf
        gmsh.model.geo.synchronize() # we have to sync before calling addPhysicalGroup
        
        surface_tag = 200
        gmsh.model.addPhysicalGroup(2,[surf, out[0][1], out[2][1], out[3][1], out[4][1], out[5][1]],surface_tag)
        gmsh.model.setPhysicalName(2,surface_tag,self.surfName)
        
        volume_tag = 300
        gmsh.model.addPhysicalGroup(3,[out[1][1]],volume_tag)
        gmsh.model.setPhysicalName(3,volume_tag,self.volName)
        
        gmsh_mesh_build("hex")
        gmsh.write(meshFileName)
        print(f"Generated {meshFileName}: {infos(gmsh.option,'hex')}")
        gmsh.finalize()
        
class Cylinder(object):
    def __init__ (self,radius,thickness,mesh_size,surfName,volName,fileFormat=4.1,verbose=False):
        """ 
            geometrical cylinder is zero centered, with radius r and length t along (Oz)
            this cylinder mesh is generated with python geo
            by default gmsh file format 4.1 is used to write the mesh text file
        """
        
        self.r = radius
        self.t = thickness
        self.msh_s = mesh_size
        self.surfName = surfName
        self.volName = volName
        self.withExtraSurf = False
        gmsh_init("cyl",fileFormat,0,verbose)

    def addEdgeSurf(self,name1,name2):
        """ optional surfaces : name1 will refer to the base surface, name2 will refer to the translated name1 surface from extrusion """
        self.n1 = name1
        self.n2 = name2
        self.withExtraSurf = True
    
    def make(self,meshFileName):
        """ write cylinder mesh file """
        
        p_origin = gmsh.model.geo.addPoint(0,0,-0.5*self.t,self.msh_s)

        start_pt = gmsh.model.geo.addPoint(self.r,0,-0.5*self.t,self.msh_s)
        interim_pt = gmsh.model.geo.addPoint(0,self.r,-0.5*self.t,self.msh_s)
        interim1_pt = gmsh.model.geo.addPoint(-self.r,0,-0.5*self.t,self.msh_s)
        interim2_pt = gmsh.model.geo.addPoint(0,-self.r,-0.5*self.t,self.msh_s)
        
        big_circle = gmsh.model.geo.addCircleArc(start_pt,p_origin,interim_pt)
        big_circle1 = gmsh.model.geo.addCircleArc(interim_pt,p_origin,interim1_pt)
        big_circle2 = gmsh.model.geo.addCircleArc(interim1_pt,p_origin,interim2_pt)
        big_circle3 = gmsh.model.geo.addCircleArc(interim2_pt,p_origin,start_pt)

        gmsh.model.geo.remove([(p_origin,0)])

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
        
        gmsh.write(meshFileName)
        print(f"Generated {meshFileName}: {infos(gmsh.option,'cyl')}")
        gmsh.finalize()

class Tube(object):
    def __init__ (self,radius1,radius2,length,mesh_size,surfName,volName,fileFormat=4.1,verbose=False):
        """
            geometrical Tube is zero centered, with inner radius radius1 and outer radius
            radius2 and length along (Oz)
            this tube mesh is generated with python geo
            by default gmsh file format 4.1 is used to write the mesh text file
        """

        self.r1 = radius1
        self.r2 = radius2
        self.l = length
        self.msh_s = mesh_size
        self.surfName = surfName
        self.volName = volName
        self.withExtraSurf = False
        gmsh_init("tube",fileFormat,0,verbose)

    def addEdgeSurf(self,name1,name2):
        """ optional surfaces : name1 will refer to the base surface, name2 will refer to the translated name1 surface from extrusion """
        self.n1 = name1
        self.n2 = name2
        self.withExtraSurf = True

    def make(self,meshFileName):
        """ write tube mesh file """

        p_origin = gmsh.model.geo.addPoint(0,0,-0.5*self.l,self.msh_s)

        start_pt = gmsh.model.geo.addPoint(self.r2,0,-0.5*self.l,self.msh_s)
        interim_pt = gmsh.model.geo.addPoint(0,self.r2,-0.5*self.l,self.msh_s)
        interim1_pt = gmsh.model.geo.addPoint(-self.r2,0,-0.5*self.l,self.msh_s)
        interim2_pt = gmsh.model.geo.addPoint(0,-self.r2,-0.5*self.l,self.msh_s)

        start_pt2 = gmsh.model.geo.addPoint(self.r1,0,-0.5*self.l,self.msh_s)
        interim3_pt = gmsh.model.geo.addPoint(0,-self.r1,-0.5*self.l,self.msh_s)
        interim4_pt = gmsh.model.geo.addPoint(-self.r1,0,-0.5*self.l,self.msh_s)
        interim5_pt = gmsh.model.geo.addPoint(0,self.r1,-0.5*self.l,self.msh_s)

        big_circle = gmsh.model.geo.addCircleArc(start_pt,p_origin,interim_pt)
        big_circle1 = gmsh.model.geo.addCircleArc(interim_pt,p_origin,interim1_pt)
        big_circle2 = gmsh.model.geo.addCircleArc(interim1_pt,p_origin,interim2_pt)
        big_circle3 = gmsh.model.geo.addCircleArc(interim2_pt,p_origin,start_pt)

        small_circle = gmsh.model.geo.addCircleArc(start_pt2,p_origin,interim3_pt)
        small_circle1 = gmsh.model.geo.addCircleArc(interim3_pt,p_origin,interim4_pt)
        small_circle2 = gmsh.model.geo.addCircleArc(interim4_pt,p_origin,interim5_pt)
        small_circle3 = gmsh.model.geo.addCircleArc(interim5_pt,p_origin,start_pt2)

        gmsh.model.geo.remove([(p_origin,0)])

        curvedLoop = gmsh.model.geo.addCurveLoop([big_circle,big_circle1,big_circle2,big_circle3]) # curvedLoop is an index (integer)
        #the following curve defines the hole
        curvedLoop2 = gmsh.model.geo.addCurveLoop([small_circle,small_circle1,small_circle2,small_circle3])

        surf = gmsh.model.geo.addPlaneSurface([curvedLoop,curvedLoop2]) # surf is an index (integer)
        out = gmsh.model.geo.extrude([(2,surf)],0,0,self.l) # 2 is the dimension of the object refered by index surf
        gmsh.model.geo.synchronize() # we have to sync before calling addPhysicalGroup

        surface_tag = 200
        gmsh.model.addPhysicalGroup(2,[surf, out[0][1], out[2][1], out[3][1], out[4][1], out[5][1],
        out[6][1], out[7][1], out[8][1], out[9][1] ],surface_tag)
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

        gmsh.write(meshFileName)
        print(f"Generated {meshFileName}: {infos(gmsh.option,'tube')}")
        gmsh.finalize()

class Ellipsoid(object):
    def __init__ (self,r1,r2,mesh_size,surfName,volName,fileFormat=4.1,verbose=False):
        """ 
            geometrical ellipsoid is zero centered, with radius r1 in the plane (Oxy) and r2 along (Oz)
            this ellipsoid mesh is generated with open cascade (occ)
            by default gmsh file format 4.1 is used to write the mesh text file
        """
        
        self.r1 = r1
        self.r2 = r2
        self.msh_s = mesh_size
        self.surfName = surfName
        self.volName = volName
        gmsh_init("ellipsoid",fileFormat,0,verbose)

    def make(self,meshFileName):
        """ write ellipsoid mesh file """
        x = 0.0
        y = 0.0
        z = 0.0
        a = 1.0
        b = 1.0
        c = self.r2/self.r1
        sph = gmsh.model.occ.addSphere(x,y,z,self.r1)
        gmsh.model.occ.dilate([(3,sph)],x,y,z,a,b,c)
        gmsh.model.occ.synchronize() # we have to sync before calling addPhysicalGroup

        surfaces = gmsh.model.getEntities(dim=2)
        surface_tag = 200
        gmsh.model.addPhysicalGroup(surfaces[0][0], [surfaces[0][1]], surface_tag) # first parameter should be 2
        gmsh.model.setPhysicalName(2,surface_tag,self.surfName)

        volumes = gmsh.model.getEntities(dim=3)
        volume_tag = 300
        gmsh.model.addPhysicalGroup(volumes[0][0], [volumes[0][1]], volume_tag) # first parameter should be 3
        gmsh.model.setPhysicalName(3,volume_tag,self.volName)

        gmsh.model.occ.synchronize() # we have to synchronize before the 'generate' call to build the mesh
        eps = 1e-3
        rescaler = 1.0 + eps
        xmin = -rescaler*self.r1
        ymin = -rescaler*self.r1
        zmin = -rescaler*self.r2
        xmax = rescaler*self.r1
        ymax = rescaler*self.r1
        zmax = rescaler*self.r2
        objects = gmsh.model.getEntitiesInBoundingBox(xmin,ymin,zmin,xmax,ymax,zmax, dim = 0)
        gmsh.model.mesh.setSize(objects,self.msh_s)
        gmsh.model.mesh.generate(dim = 3)
        gmsh.write(meshFileName)

        print(f"Generated {meshFileName}: {infos(gmsh.option,'ellipsoid')}")
        gmsh.finalize()

# Cuboid class
# make builds a mesh of a cuboid format 2.2. Each elementary cube split in 6 tetrahedrons

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
        """ write mesh file in gmsh 2.2 text format """
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
        print(f"Generated {meshFileName}: {len(self.pts)} nodes, {len(self.Tet)} tetrahedra")

# Ico class: an icosahedron mesh from scratch
class Ico(object):
    def __init__ (self):
        """
            builds a mesh of an icosahedron in text format 2.2, each tetrahedron is built with a outer face of the icosahedron linked to the barycenter.
            see https://fr.wikipedia.org/wiki/Icosa%C3%A8dre for more informations. Mesh file is written in the output text file without gmsh module.
        """
        self.phi = (1+sqrt(5))/2
        self.pts = np.array([[self.phi,1,0], [self.phi,-1,0], [-self.phi,1,0], [-self.phi,-1,0],[1,0,self.phi], [1,0,-self.phi], [-1,0,self.phi], [-1,0,-self.phi],[0,self.phi,1], [0,self.phi,-1], [0,-self.phi,1], [0,-self.phi,-1], [0,0,0]] )

        self.f_idx = [ [0,8,4], [0,4,1], [0,1,5], [0,5,9], [0,9,8], [8,9,2], [8,2,6], [8,6,4], [6,10,4], [6,3,10], [6,2,3], [10,3,11], [10,11,1], [10,1,4], [1,11,5], [11,3,7], [11,7,5], [3,2,7], [7,9,2],  [5,7,9] ]

    def make(self,meshFileName,volRegionName,surfRegionName):
        """ write mesh file in gmsh 2.2 text format """
        meshFile = open(meshFileName,'w')
        meshFile.write("$MeshFormat\n2.2\t0\t8\n$EndMeshFormat\n")
        meshFile.write("$PhysicalNames\n2\n")
        surfRegionTag = 200
        volRegionTag = 300
        meshFile.write("2 " + str(surfRegionTag) + ' \"' + surfRegionName + '\"\n')
        meshFile.write("3 " + str(volRegionTag) + ' \"' + volRegionName + '\"\n')
        meshFile.write("$EndPhysicalNames\n")
        meshFile.write("$Nodes\n" + str(len(self.pts))+"\n")
        for i in range(0,len(self.pts)):
            meshFile.write( str(i+1)+"\t"+str(self.pts[i][0])+"\t"+str(self.pts[i][1])+"\t"+str(self.pts[i][2]) + "\n" )
        meshFile.write("$EndNodes\n$Elements\n")
        meshFile.write(str(2*len(self.f_idx))+"\n")
        for i in range(0,len(self.f_idx)):
            meshFile.write( str(i+1)+"\t4\t2\t" + str(volRegionTag) + "\t1\t13\t"+str(1+self.f_idx[i][0])+"\t"+str(1+self.f_idx[i][1])+"\t"+str(1+self.f_idx[i][2]) + "\n" )
        for i in range(0,len(self.f_idx)):
            meshFile.write( str(len(self.f_idx)+i+1)+"\t2\t2\t" + str(surfRegionTag) + "\t1\t"+str(1+self.f_idx[i][0])+"\t"+str(1+self.f_idx[i][1])+"\t"+str(1+self.f_idx[i][2]) + "\n" )

        meshFile.write("$EndElements\n")
        meshFile.close()
        print(f"Generated {meshFileName}: 13 nodes, 20 tetrahedra, 20 triangles")

