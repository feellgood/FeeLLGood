rev_CMK is a command line tool to apply the Cuthill McKee algorithm (CMK) to a mesh file, format gmsh 2.2
revCMK requires BOOST. 
To build the executable, type from a terminal
cmake .
make

 you may eventually install the executable typing
sudo make install

Then to apply CMK algorithm to a .msh file, type from command line:
rev_CMK my_file.msh

The input file 'my_file.msh' is not modified.
The application rev_CMK will create a new file with an extra extension .r_cmk, here it would have been 'my_file.msh.r_cmk'
If ever some extra informations are stored in the input file, they are copied to the output file. It may be convenient for some metadata, or comments. 

