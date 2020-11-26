# to install vtk for python type the following command:
# python -m pip install --upgrade --user vtk

import vtk
import sys
import os
import numpy as np

class mesh(object):
    def __init__ (self,fileName):
        self.Nodes = []
        self.Tet = []
        mshFile = open(fileName,'r')
        count = 0
        current_line =''
        while current_line.strip() != '$MeshFormat':
            current_line = mshFile.readline()

        mshFormat = mshFile.readline()
        if mshFormat.strip() != "2.2\t0\t8\n":
            print("mesh format = 2.2")
        else:
            print("cannot read " + fileName + "file, wrong format, 2.2 required.")
            sys.exit(1)

        while current_line.strip() != '$Nodes':
            current_line = mshFile.readline()
        current_line = mshFile.readline()
        
        while current_line.strip() != '$EndNodes':
            current_line = mshFile.readline()
            ls = current_line.strip().split("\t")
            if len(ls) == 4:
                self.Nodes.append( [float(ls[1]),float(ls[2]),float(ls[3]) ])
        current_line = mshFile.readline()
        current_line = mshFile.readline()
        
        while current_line.strip() != '$EndElements':
            current_line = mshFile.readline()
            ls = current_line.strip().split("\t")
            if len(ls) == 9: # we only consider tetrahedrons, assuming there is no other objects than triangles and tetrahedrons in the Elements list
                self.Tet.append( [int(ls[5])-1,int(ls[6])-1,int(ls[7])-1,int(ls[8])-1 ])
        mshFile.close()
        print("mesh read from file " + fileName)

    def infos(self):
        print("nb Nodes: ",len(self.Nodes),"\nnb Tetrahedrons: ",len(self.Tet))

class mag(object):
    def __init__(self,fileName):
        self.Mag = vtk.vtkFloatArray()
        self.Mag.SetNumberOfComponents(3)
        magFile = open(fileName,'r')
        lines = magFile.readlines()        
        nb_values =len(lines) # first .sol line is always "#t = xxx"
        self.Mag.SetNumberOfTuples(nb_values)
        
        for i in range(1,nb_values):
            ls = lines[i].strip().split('\t')
            self.Mag.SetTuple3(i-1,float(ls[4]),float(ls[5]),float(ls[6]))
        magFile.close()

def get_params():
    import argparse
    description = 'convert .sol to .vtk using gmsh mesh.'
    epilogue = '''
 convert feellgood output .sol files to .vtk file using mesh in gmsh version 2.2 format 
   '''
    parser = argparse.ArgumentParser(description=description, epilog=epilogue,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('mshFileName', help='truc.msh')
    parser.add_argument('solFileName', help='truc.sol')
    args = parser.parse_args()
    return [args.mshFileName,args.solFileName]

def main():
    print("This is convert2vtk python script, using vtk version " + vtk.vtkVersion.GetVTKVersion())
    FileNames = get_params()
    path, ext = os.path.splitext(FileNames[0])
    ext = ext.lower()
    if ext == '.msh':
        print("input mesh file is ", FileNames[0])
    
    print("input sol file name =", FileNames[1])
    
    msh = mesh(FileNames[0])
    msh.infos()

    data = mag(FileNames[1])

    points = vtk.vtkPoints()
    points.SetNumberOfPoints(len(msh.Nodes))

    for i in range(0,len(msh.Nodes)):
        points.InsertPoint(i,msh.Nodes[i])
    ugrid = vtk.vtkUnstructuredGrid()
    #ugrid.Allocate(len(msh.Tet))
    ugrid.SetPoints(points)
    for i in range(0,len(msh.Tet)):
        ugrid.InsertNextCell(vtk.VTK_TETRA,4,msh.Tet[i])

    ugrid.GetPointData().SetVectors(data.Mag)
    ugrid.Modified()# usefull ?

    writer = vtk.vtkXMLUnstructuredGridWriter()# to replace by vtkGenericDataObjectWriter ? or vtkDataSetWriter ? eventually non XML ?
    writer.SetFileName('toto.vtk')
    writer.SetInputData(ugrid)
    writer.Update()# usefull ?
    writer.SetDataModeToAscii() #only if XML
    writer.Write()
    
if __name__ == "__main__":
    main()
