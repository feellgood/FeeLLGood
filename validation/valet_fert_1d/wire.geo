// Wire of square section a x a and length L along z, for the 1D validation of the spin injection
// through a polarizer stack (valet_fert.py). Lengths in nm (length_unit: 1e-9).
SetFactory("OpenCASCADE");
L = 300;
a = 10;
h = 1;
Box(1) = {-a/2, -a/2, 0, a, a, L};
Mesh.CharacteristicLengthMin = h;
Mesh.CharacteristicLengthMax = h;
Mesh.MshFileVersion = 4.1;
// OpenCASCADE box: surfaces 1, 2 at x min, max; 3, 4 at y min, max; 5, 6 at z min, max
Physical Volume("wire") = {1};
Physical Surface("inlet") = {5};
Physical Surface("outlet") = {6};
Physical Surface("side") = {1, 2, 3, 4};
