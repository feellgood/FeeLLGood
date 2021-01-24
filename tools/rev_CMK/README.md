`rev_CMK` is a command line tool that attempts to reduce the graph
bandwidth of a mesh by reordering its nodes. It implements the reverse
Cuthill–McKee algorithm. Its input file should be a mesh in the GMSH 2.2
format. `rev_CMK` requires the Boost library.

To build the executable, type in a terminal

```shell
cmake .
make
```

you may install it by typing

```shell
sudo make install
```

Then, to process a .msh file, type:

```shell
rev_CMK my_file.msh
```

The input file `my_file.msh` is not modified. The output will be written
to a new file named with the extra extension `.r_cmk`, in this case
`my_file.msh.r_cmk`.

If the input file contains extra information, such as metadata or
comments, that information is copied to the output.
