rev_CMK is a command line tool that applies the Cuthill McKee algorithm
(CMK) to a mesh file, format gmsh 2.2. revCMK requires Boost. 

To build the executable, type from a terminal

    cmake .
    make

you may install the executable by typing

    sudo make install

Then, to apply the CMK algorithm to a .msh file, type from the command
line:

    rev_CMK my_file.msh

The input file 'my_file.msh' is not modified.

The application rev_CMK will create a new file with an extra extension
.r_cmk. Here it would have been 'my_file.msh.r_cmk'

If any extra information is stored in the input file, it is copied to
the output file. It may be convenient for some metadata, or comments.
