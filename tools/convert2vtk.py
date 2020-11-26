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
        self.Mag = []
        magFile = open(fileName,'r')
        current_line = magFile.readline()
        while True:
            current_line = magFile.readline()
            ls = current_line.strip().split('\t')
            if len(ls) >= 7:
                self.Mag.append([float(ls[4]),float(ls[5]),float(ls[6])])
            if not current_line:
                break
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
    ugrid.Modified()# usefull ?

    vects = vtk.vtkFloatArray()
    vects.SetNumberOfComponents(3)
    vects.SetNumberOfTuples(len(data.Mag))
    for i in range(0,len(data.Mag)):
        vects.SetTuple3(i,data.Mag[i][0],data.Mag[i][1],data.Mag[i][2])

    ugrid.GetPointData().SetVectors(vects)
    ugrid.Modified()# usefull ?

    writer = vtk.vtkXMLUnstructuredGridWriter()# to replace by vtkGenericDataObjectWriter ? or vtkDataSetWriter ? eventually non XML ?
    writer.SetFileName('toto.vtk')
    writer.SetInputData(ugrid)
    writer.Update()# usefull ?
    writer.SetDataModeToAscii() #only if XML
    writer.Write()
    
    arrow = vtk.vtkArrowSource()
    arrow.SetTipResolution(16)
    arrow.SetTipLength(0.3)
    arrow.SetTipRadius(0.1)

    glyph = vtk.vtkGlyph3D()
    glyph.SetInputData(ugrid)
    glyph.SetSourceConnection(arrow.GetOutputPort())
    glyph.SetVectorModeToUseVector()
    glyph.SetColorModeToColorByVector()
    glyph.OrientOn()
    glyph.Update()

    mbds = vtk.vtkMultiBlockDataSet()
    mbds.SetNumberOfBlocks(2)
    mbds.SetBlock(0, ugrid)
    mbds.SetBlock(1, glyph.GetOutput())

    mapper = vtk.vtkCompositePolyDataMapper2()
    mapper.SetInputDataObject(mbds)
    cdsa = vtk.vtkCompositeDataDisplayAttributes()
    mapper.SetCompositeDataDisplayAttributes(cdsa)


#    ugridMapper = vtk.vtkPolyDataMapper() #vtk.vtkDataSetMapper()
#    ugridMapper.SetInputConnection(glyph.GetOutputPort()) #SetInputData(ugrid)
#    ugridMapper.ScalarVisibilityOn()
#    ugridMapper.SetScalarModeToUsePointData()
#    ugridMapper.SetScalarRange(ugrid.GetPointData().GetVectors().GetRange(-1))

    ugridActor = vtk.vtkActor()
    ugridActor.SetMapper(mapper)#(ugridMapper)
    colors = vtk.vtkNamedColors()
    ugridActor.GetProperty().SetColor(colors.GetColor3d("Peacock"))
    ugridActor.GetProperty().EdgeVisibilityOn()

    renderer = vtk.vtkRenderer()
    renderer.AddActor(ugridActor)
    
    renderer.SetBackground(colors.GetColor3d("Beige"))

    renderer.ResetCamera()
    renderer.GetActiveCamera().Elevation(60.0)
    renderer.GetActiveCamera().Azimuth(30.0)
    renderer.GetActiveCamera().Dolly(1.2)

    

    renWin = vtk.vtkRenderWindow()
    renWin.AddRenderer(renderer)
    iren = vtk.vtkRenderWindowInteractor()
    iren.SetRenderWindow(renWin)
    renWin.SetSize(640, 480)
    renWin.SetWindowName(FileNames[0])
    # Interact with the data.
    renWin.Render()

    iren.Start()

if __name__ == "__main__":
    main()
