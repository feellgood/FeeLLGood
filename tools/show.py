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
    print("This is show python script, using vtk version " + vtk.vtkVersion.GetVTKVersion())
    fileName = get_params()
    path, ext = os.path.splitext(fileName)
    ext = ext.lower()
    if ext == '.vtk':
        print("input file is ", fileName)
    else:
        print(ext, " unsupported.")
        sys.exit(1)

    reader = vtk.vtkXMLUnstructuredGridReader()
    reader.SetFileName(fileName)
    reader.Update()

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
