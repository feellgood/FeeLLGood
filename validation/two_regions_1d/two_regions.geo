// Wire of square section a x a along z: non-magnetic spacer N (0 < z < tN) and free layer F2
// (tN < z < L), for the validation of the continuity of s between regions of different densities of
// states (two_regions.py). Lengths in nm (length_unit: 1e-9).
SetFactory("OpenCASCADE");
L = 300;
tN = 40;
a = 4;
h = 1;
Box(1) = {-a/2, -a/2, 0, a, a, tN};
Box(2) = {-a/2, -a/2, tN, a, a, L - tN};
BooleanFragments{ Volume{1}; Delete; }{ Volume{2}; Delete; }
Mesh.CharacteristicLengthMin = h;
Mesh.CharacteristicLengthMax = h;
Mesh.MshFileVersion = 4.1;
e = 1e-3;
Physical Volume("N") = Volume In BoundingBox{-a, -a, -e, a, a, tN + e};
Physical Volume("F2") = Volume In BoundingBox{-a, -a, tN - e, a, a, L + e};
Physical Surface("inlet") = Surface In BoundingBox{-a, -a, -e, a, a, e};
Physical Surface("outlet") = Surface In BoundingBox{-a, -a, L - e, a, a, L + e};
