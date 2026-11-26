# feeLLGood tools in C++/Qt

C++ versions of the Python tools of the parent directory, built with Qt 5.
They take the same arguments and write the same files, except where noted
below.

* **fg-3D**: 3D viewer for a gmsh mesh and a `.sol` file: surface,
  magnetization arrows, surface charges, export to `.vtu`
* **fg-gensol**: generate an initial magnetization file `sol.in` from a
  mesh and a function per volume region
* **fg-merge-sol**: combine a mesh and a `.sol` file into a table of the
  physical quantities as a function of position
* **fg-plot**: quiver plot of a `.sol` file, in 3D or in a 2D subspace
* **fg-sol2hdf5**: archive the `.sol`, `.evol`, settings and mesh files of
  a simulation in a HDF5 file
* **fg-sol2vtk**: convert a `.sol` file and its mesh to ParaView's `.vtu`
  format

Meshes are read with the gmsh API, as feeLLGood does. HDF5 files are
interchangeable with the ones of the Python tools.

## Dependencies

* Qt 5 (Core, Gui, Widgets, Qml), the version VTK's `GUISupportQt` module
  is built with
* gmsh (library and `gmsh.h`), as for feeLLGood
* VTK 9 with its Qt support, for fg-3D, fg-plot and fg-sol2vtk
* HDF5 with its C++ API, for fg-sol2hdf5
* yaml-cpp, as for feeLLGood, for fg-sol2hdf5

A tool whose optional dependency is missing is not built. On Debian or
Ubuntu:

```shell
sudo apt install qtbase5-dev qtdeclarative5-dev libgmsh-dev libvtk9-dev \
    libvtk9-qt-dev libhdf5-dev libyaml-cpp-dev
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

* **fg-gensol** reads the magnetization functions from `init_cfg.js`
  instead of `init_cfg.py`. It is written on the first run and maps each
  region name to a JavaScript function:

  ```js
  ({
      'region': function (x, y, z) { return [cos(x), sin(x), 0]; },
  })
  ```

  The functions of the JavaScript `Math` object can be used without the
  `Math.` prefix. If `init_cfg.js` has an error or lacks a region,
  fg-gensol prints the problem, or the lines to add, and leaves the file
  unchanged.
* **fg-plot** evaluates the `--filter` predicate as JavaScript.
  Python's `and`, `or`, `not`, `abs`, `pi` and the functions of the `math`
  module work, but chained comparisons do not: write `0 < x && x < 1`,
  not `0 < x < 1`. The plots cannot be zoomed or panned as with
  matplotlib, but the 3D view can be rotated.
* **fg-sol2vtk**, **fg-plot** and **fg-merge-sol** check that the `.sol`
  file has as many data lines as the mesh has nodes.
* **fg-sol2hdf5** checks that a `.evol` file exists before writing
  anything.
* **fg-3D** checks that the `.sol` file matches the mesh before
  exporting to `.vtu`.
