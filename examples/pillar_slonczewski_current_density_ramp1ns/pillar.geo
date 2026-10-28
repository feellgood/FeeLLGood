// Free layer F2 of the example pillar_slonczewski_current_density_ramp1ns: disk of radius 20 nm and
// thickness 4 nm, axis z, z in [-2, 2] nm, mesh size 2 nm (the s-d length l_sd = 1 nm is not
// resolved: see the README).
// Built-in geometry kernel of gmsh, as meshMaker.Cylinder: a disk of four circle arcs, extruded
// along z.
// Surfaces: pillar_bottom (z = -2 nm, current inlet, polarizer stack), pillar_top (z = +2 nm,
// outlet), pillar_side; volume: F2. Lengths in nm (settings: mesh.length_unit = 1e-9).
//
//   gmsh -3 pillar.geo -o pillar.msh
R = 20;
t = 4;
h = 2;

Point(1) = {0, 0, -t/2, h};     // centre of the arcs, deleted below
Point(2) = {R, 0, -t/2, h};
Point(3) = {0, R, -t/2, h};
Point(4) = {-R, 0, -t/2, h};
Point(5) = {0, -R, -t/2, h};
Circle(1) = {2, 1, 3};
Circle(2) = {3, 1, 4};
Circle(3) = {4, 1, 5};
Circle(4) = {5, 1, 2};
Delete { Point{1}; }            // not part of the geometry (no isolated node in the mesh)
Curve Loop(1) = {1, 2, 3, 4};
Plane Surface(1) = {1};

// out[0]: top surface, out[1]: volume, out[2..5]: lateral surfaces
out[] = Extrude {0, 0, t} { Surface{1}; };

Physical Volume("F2") = {out[1]};
Physical Surface("pillar_bottom") = {1};
Physical Surface("pillar_top") = {out[0]};
Physical Surface("pillar_side") = {out[2], out[3], out[4], out[5]};

Mesh.MshFileVersion = 4.1;
