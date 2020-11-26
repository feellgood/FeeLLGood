# to install vtk for python type the following command:
# python -m pip install --upgrade --user vtk

import vtk
import sys
import os

def get_params():
    import argparse
    description = 'make a rendering of a .vtk file'
    epilogue = '''
 render a .vtk file, assuming blablabla
   '''
    parser = argparse.ArgumentParser(description=description, epilog=epilogue,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('vtkFileName', help='truc.vtk')
    args = parser.parse_args()
    return args.vtkFileName

def main():
    fileName = get_params()
    print("This is show python script, using vtk version " + vtk.vtkVersion.GetVTKVersion()," rendering ", fileName)
    path, ext = os.path.splitext(fileName)
    ext = ext.lower()
    if ext != '.vtk':
        print(ext, " unsupported.")
        sys.exit(1)

    reader = vtk.vtkXMLUnstructuredGridReader()
    reader.SetFileName(fileName)
    reader.Update()
    output = reader.GetOutput()
#    scalar_range = output.GetScalarRange()

#    ugrid = vtk.vtkUnstructuredGrid()

#    it = reader.GetOutput().NewCellIterator()
#    it.InitTraversal()
#    while not it.IsDoneWithTraversal():
#        cell = vtk.vtkGenericCell()
#        it.GetCell(cell)
        
#        cellName = vtk.vtkCellTypes.GetClassNameFromTypeId(cell.GetCellType())
#        print(cellName, "NumberOfPoints:", cell.GetNumberOfPoints(), "CellDimension:", cell.GetCellDimension())
#        print("cell_pts = ",cell.GetPoints())
#        if cell.GetNumberOfPoints() == 4:
#            ugrid.InsertNextCell(vtk.VTK_TETRA,4,cell.GetPoints()) # ouch le cast sur tetraèdre
#        it.GoToNextCell()

    arrow = vtk.vtkArrowSource()
    arrow.SetTipResolution(16)
    arrow.SetTipLength(0.3)
    arrow.SetTipRadius(0.1)
    arrow.SetShaftResolution(16)
    arrow.SetShaftRadius(0.03)

    glyph = vtk.vtkGlyph3D()
    glyph.SetInputData(reader.GetOutput())#(ugrid)
    glyph.SetSourceConnection(arrow.GetOutputPort())
    glyph.SetVectorModeToUseVector() # to use the vector input data
    glyph.SetColorModeToColorByVector() #SetColorModeToColorByScalar()
    glyph.SetScaleModeToDataScalingOff()
    glyph.OrientOn()
    glyph.Update()

#    mbds = vtk.vtkMultiBlockDataSet()
#    mbds.SetNumberOfBlocks(2)
#    mbds.SetBlock(0, ugrid)
#    mbds.SetBlock(1, glyph.GetOutput())

#    mapper = vtk.vtkDataSetMapper()
#    mapper.SetInputDataObject(output)#(mbds)


    mapper = vtk.vtkPolyDataMapper() #vtk.vtkDataSetMapper()
    mapper.SetInputConnection(glyph.GetOutputPort()) #SetInputData(ugrid)
    mapper.ScalarVisibilityOn()
    mapper.SetScalarModeToUsePointData()
#    mapper.SetScalarRange(ugrid.GetPointData().GetVectors().GetRange(-1))

    actor = vtk.vtkActor()
    actor.SetMapper(mapper)#(ugridMapper)
    colors = vtk.vtkNamedColors()
    actor.GetProperty().SetColor(colors.GetColor3d("Peacock"))
    actor.GetProperty().EdgeVisibilityOn()

    renderer = vtk.vtkRenderer()
    renderer.AddActor(actor)
    
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
    renWin.SetWindowName(fileName)
    # Interact with the data.
    renWin.Render()

    iren.Start()

if __name__ == "__main__":
    main()
