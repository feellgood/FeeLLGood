Those two python scripts use VTK, convert2vtk.py is a command line tool to convert a feellgood .sol output file to a .vtu file. You need to provide the mesh file.
To install VTK for python, type in a terminal:

python -m pip install --upgrade --user vtk

How to use the script:
Example, from /tools subdirectory:

./convert2vtk mymesh.msh mysol.sol
The script will generate a file named mysol.vtu, using XML format and binary for the datas.
The file mysol.vtu can be directly opened with paraview, it is also readable using XML Unstructured Grid Reader in VTK.

Some possible use of .vtu files in paraview:
to make some iso-surface of a component of the magnetization field, use the filter 'ExtractComponent' and build the iso-surface using the filter 'Contour' on the extracted component. 

if you have installed feellgood using 'sudo make install' 
and typed
sudo chmod +x /usr/local/bin/convert2vtk

Then you can launch convert2vtk from any directory
