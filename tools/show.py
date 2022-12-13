#!/usr/bin/env python3
# to install vtk for python type the following command:
# python -m pip install --upgrade --user vtk

import vtk
import sys
import os

def get_params():
    import argparse
    description = 'make a rendering of a .vtu file'
    epilogue = '''
 render a .vtu file assuming the input file is a XML unstructured grid file
   '''
    parser = argparse.ArgumentParser(description=description, epilog=epilogue,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('fileName', help='render a .vtu file')
    
    parser.add_argument('-P', '--percent',
    default=100,
    help='Percentage of the field vectors rendered'
    )
    args = parser.parse_args()
    
    return args

def main():
    my_args = get_params()
    fileName = my_args.fileName
    perc = my_args.percent
    #print( "fileName = ", fileName , "percentage = ", perc  )
    print("This is show python script, using VTK version " + vtk.vtkVersion.GetVTKVersion()," rendering ", fileName)
    path, ext = os.path.splitext(fileName)
    ext = ext.lower()
    if ext != '.vtu':
        print(ext, " unsupported.")
        sys.exit(1)

    reader = vtk.vtkXMLUnstructuredGridReader()
    reader.SetFileName(fileName)
    reader.Update()
    output = reader.GetOutput()
    #print(output)

    if(float(perc) < 100):
        print("let's subsample!") # we could use a vtkTransformFilter to subsample
#    contour = vtk.vtkContourFilter()
#    contour.SetInputData(output)
#    contour.ComputeNormalsOn()
#    contour.SetValue(0,0.0)
#    contour.Update()

#    isoMapper = vtk.vtkPolyDataMapper()
#    isoMapper.SetInputData(contour.GetOutput())
#    isoMapper.ScalarVisibilityOn()

# Take the isosurface data and create geometry
#actorIso = vtk.vtkLODActor()
#actorIso.SetNumberOfCloudPoints( 1000000 )
#actorIso.SetMapper( isoMapper )
#actorIso.GetProperty().SetColor( 1, 1, 1 ) 

    arrow = vtk.vtkArrowSource()
    arrow.SetTipResolution(16)
    arrow.SetTipLength(0.3)
    arrow.SetTipRadius(0.1)
    arrow.SetShaftResolution(16)
    arrow.SetShaftRadius(0.03)

    glyph = vtk.vtkGlyph3D()
    glyph.SetInputData(output)  #(reader.GetOutput())#(ugrid)
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
    mapper.SetInputConnection(glyph.GetOutputPort())#(contour.GetOutputPort())
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
    renWin.SetSize(1000, 1000)
    renWin.SetWindowName(fileName)
    # Interact with the data.
    renWin.Render()

    iren.Start()

if __name__ == "__main__":
    main()
