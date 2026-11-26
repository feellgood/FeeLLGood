# feeLLGood tools in C++/Qt

C++ versions of the Python tools of the parent directory, built with Qt 5.
They take the same arguments and write the same files, except where noted
below.

* **fg-3D**: 3D viewer for a gmsh mesh and a `.sol` file: surface,
  magnetization arrows, surface charges, export to `.vtu`
* **fg-merge-sol**: combine a mesh and a `.sol` file into a table of the
  physical quantities as a function of position
* **fg-sol2vtk**: convert a `.sol` file and its mesh to ParaView's `.vtu`
  format

Meshes are read with the gmsh API, as feeLLGood does.

## Dependencies

* Qt 5 (Core, Gui, Widgets), the version VTK's `GUISupportQt` module
  is built with
* gmsh (library and `gmsh.h`), as for feeLLGood
* VTK 9 with its Qt support, for fg-3D and fg-sol2vtk

A tool whose optional dependency is missing is not built. On Debian or
Ubuntu:

```shell
sudo apt install qtbase5-dev libgmsh-dev libvtk9-dev libvtk9-qt-dev
```

## Build

```shell
cmake -S . -B build
cmake --build build
```

`cmake --install build` installs the tools in the `bin` directory of the
install prefix. They have the names of the Python tools, which feeLLGood's
`make install` also installs: use another prefix to keep both.

## Differences with the Python tools

* **fg-sol2vtk** and **fg-merge-sol** check that the `.sol`
  file has as many data lines as the mesh has nodes.
* **fg-3D** checks that the `.sol` file matches the mesh before
  exporting to `.vtu`.
