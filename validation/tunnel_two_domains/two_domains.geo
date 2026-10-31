// Slab F2 of size Lx x Ly x t (x along the domain wall normal, z along the current), for the
// validation of the tunnel barrier with two magnetic domains (two_domains.py). Lengths in nm
// (length_unit: 1e-9).
SetFactory("OpenCASCADE");
Lx = 100;
Ly = 2;
t = 20;
h = 1;
Box(1) = {-Lx/2, -Ly/2, 0, Lx, Ly, t};
Mesh.CharacteristicLengthMin = h;
Mesh.CharacteristicLengthMax = h;
Mesh.MshFileVersion = 4.1;
// OpenCASCADE box: surfaces 1, 2 at x min, max; 3, 4 at y min, max; 5, 6 at z min, max
Physical Volume("F2") = {1};
Physical Surface("inlet") = {5};
Physical Surface("outlet") = {6};
Physical Surface("side") = {1, 2, 3, 4};
