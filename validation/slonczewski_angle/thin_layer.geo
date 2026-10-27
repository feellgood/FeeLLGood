// Thin free layer F2, of square section a x a and thickness t along z, for the angular validation
// of the spin injection through a polarizer stack (slonczewski_angle.py). Lengths in nm
// (length_unit: 1e-9). The mesh size h resolves the s-d length l_sd = 1 nm.
SetFactory("OpenCASCADE");
t = 3;
a = 2;
h = 0.2;
Box(1) = {-a/2, -a/2, 0, a, a, t};
Mesh.CharacteristicLengthMin = h;
Mesh.CharacteristicLengthMax = h;
Mesh.MshFileVersion = 4.1;
// OpenCASCADE box: surfaces 1, 2 at x min, max; 3, 4 at y min, max; 5, 6 at z min, max
Physical Volume("F2") = {1};
Physical Surface("inlet") = {5};
Physical Surface("outlet") = {6};
Physical Surface("side") = {1, 2, 3, 4};
