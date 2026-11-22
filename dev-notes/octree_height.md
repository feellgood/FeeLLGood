# Height of the FMM octree adapted to the particles

## Context

The demagnetizing field is computed in `src/fmm_demag.h` with the fast multipole method (FMM) of
ScalFMM 1.5 (rotation kernel, spherical harmonics truncated at $P = 9$). The particles are the
magnetic nodes (targets) and the Gauss points of the magnetic tetrahedrons and triangles (sources,
carrying the volume and surface magnetic charges), so that the number of particles is
$N = N_\text{nodes} + 5 N_\text{tet} + 3 N_\text{tri}$. Their positions are scaled into the
cube $[-1, 1]^3$, inside a box of width 2.01 centered on the origin.

ScalFMM sorts the particles in an octree of height $h$: the leaves are the cells of depth
$h-1$, there are $2^{h-1}$ of them per direction. Only the non-empty cells are created. The
potential at a target is the sum of

- a near field (P2P), computed directly from the particles of its leaf and of the neighbouring
  leaves;
- a far field, computed from multipole and local expansions (P2M, M2M, M2L, L2L, L2P) at all the
  levels of the octree.

The cost of the FMM depends strongly on the height. With a fixed height, such as the
`NbLevels = 6` used before commit 2fdda7690 ("Adapt the FMM octree height to the particles"),
the leaves of depth 5 of usual meshes (a few thousand nodes, about $10^5$ particles) hold a few
particles only: the far field of these levels, and especially the M2L transfers (up to 189
interactions per cell, each one costing $O(P^3)$ with the rotation kernel), costs much more than
the direct computation it replaces. The height is therefore adapted to the particles.

## Choice of the height

The height is chosen by `src/octree_height.h`:

> the height is the largest $h \in [3, 10]$ such that the non-empty leaves of the octree of
> height $h$ contain on average at least 150 particles.

The non-empty leaves are counted without building any octree:

1. the Morton key of each particle is computed at depth 9 (the depth of the leaves of the
   deepest octree allowed, $2^9$ cells per direction): the bits of the three cell indices are
   interleaved, so that the key of the cell containing the particle at depth $h-1$ is the key
   shifted right by $3(10-h)$ bits;
2. the keys are sorted: the particles of a same cell, at any depth, are then contiguous;
3. for $h = 10, 9, \dots, 4$, the number $L_h$ of distinct shifted keys is the number of
   non-empty leaves; the first $h$ such that $N \ge 150 L_h$ is returned, $h = 3$ otherwise.

Only the non-empty leaves are counted, so that the height adapts to the shape of the mesh: a thin
film or a long wire fills a small fraction of the cubic box, and gets a deeper octree than a
sphere with the same number of particles. The cost is a sort of the particles, done once when the
octree is built (17 ms for $1.1 \times 10^5$ particles on one thread).

The number of sub-levels of the octree (`SizeSubLevels`, a parameter of the storage of the cells
in ScalFMM, which must not exceed $h-1$) is 3, reduced to $h-1$ if needed.

### Calibration of the threshold

The threshold of 150 particles per leaf was calibrated with the same kernel (rotation, $P = 9$)
and 8 threads, on a family of meshes of a cylinder of diameter 20 nm and length 20 nm, described
by the following gmsh file: the base disc is meshed with the characteristic size `lc`, then
extruded along $z$ in `Nz` layers.

```
lc = 2;  L = 20;  D = 20;  Nz = L/4;
R = D/2;  Lo2 = L/2;
Point(1) = { 0,  0, -Lo2, lc};
Point(2) = { R,  0, -Lo2, lc};  Point(3) = { 0,  R, -Lo2, lc};
Point(4) = {-R,  0, -Lo2, lc};  Point(5) = { 0, -R, -Lo2, lc};
Circle(1) = {2, 1, 3};  Circle(2) = {3, 1, 4};
Circle(3) = {4, 1, 5};  Circle(4) = {5, 1, 2};
Line Loop(10) = {1, 2, 3, 4};
Plane Surface(1) = {10};
out0[] = Extrude {0, 0, L} { Surface{1}; Layers{Nz}; };
Physical Volume("magnet", 300) = {out0[1]};
Physical Surface("lateral", 200) = {out0[2], out0[3], out0[4], out0[5]};
Physical Surface("inlet", 201) = {1};
Physical Surface("outlet", 202) = {out0[0]};
```

The finer meshes are obtained by halving `lc` and doubling `Nz` (`Nz` = L / (2 `lc`)), then
`gmsh -3 cyl.geo`. These measurements were made with a version of the FMM interface that uses 4
Gauss points per triangle instead of 3, hence
$N = N_\text{nodes} + 5 N_\text{tet} + 4 N_\text{tri}$ in the table below. The last two columns
give, for comparison, the time with a fixed deeper octree.

| `lc` (nm) | `Nz` | nodes | tetrahedrons | triangles | particles | adapted height | FMM time | fixed height | FMM time |
|-----------|------|-------|--------------|-----------|-----------|----------------|----------|--------------|----------|
| 2         | 5    | 738   | 3180         | 744       | 19614     | 3              | 7 ms     | 8            | 260 ms   |
| 1         | 10   | 4653  | 23400        | 2840      | 133013    | 4              | 96 ms    | 7            | 1.4 s    |
| 0.5       | 20   | 33516 | 183720       | 11244     | 997092    | 5              | 665 ms   | 7            | 7.9 s    |
| 0.25      | 40   |       |              |           | ~7.5e6    | 6              | 6.1 s    | 7            | 26 s     |

(8 threads, time of one FMM call. The mesh counts were checked by regenerating the meshes with
gmsh 4.12.1, except the finest one.)

The fastest octree always has about 300 particles per non-empty leaf; one level deeper (about 35
particles per leaf) is 4 to 5 times slower, one level less (about 2000 particles per leaf) is 10
to 25 % slower. Since the number of particles per leaf is divided by about 8 from one level to
the next, the threshold 150 selects a height with 150 to 1200 particles per leaf, which contains
the optimum and excludes the costly extra level.

## Implementation in feeLLGood

- `src/octree_height.h` (no dependency on ScalFMM): `scal_fmm::mortonKey` and
  `scal_fmm::octreeHeight`, with the constants `MinLevels = 3`, `MaxLevels = 10` and
  `MinParticlesPerLeaf = 150`.
- `src/fmm_demag.h`: the class `scal_fmm::fmm` has a member `nbLevels`, computed from the
  normalized positions of the particles (`fmm::particlePositions`, same order and same scaling as
  the insertion in the tree). The tree and the kernel are built with this height in the
  initializer list; `nbLevels`, and the normalization coefficient `norm` that
  `particlePositions` needs, are declared before them, so that they are initialized first.
- `src/main.cpp`: in verbose mode, the height is printed with the insertion of the particles:
  `Magnetostatics: particles inserted, using 8 threads, octree of 4 levels, in 28.21 ms`.
- `unit-tests/ut_octree_height.cpp`: Morton keys of the 8 cells of depth 1, uniform grid of
  $64^3$ particles filling the box (height 4), the same grid in one cell of depth 4 (height 8:
  the empty leaves are ignored), few particles (height 3), all particles at the same point
  (height 10).

## Benchmark meshes

The measurements below use three ellipsoids of revolution generated by the class `Ellipsoid` of
`python-modules/meshMaker.py`: semi-axis `r1` in the plane $(Oxy)$, `r2` along $(Oz)$,
characteristic mesh size `mesh_size`, all in nm (`length_unit: 1e-9`). For instance, E1 is
generated by

```python
from feellgood.meshMaker import Ellipsoid

mesh = Ellipsoid(r1 = 4, r2 = 8, mesh_size = 0.5,
                 surfName = "ellipsoid_surface", volName = "ellipsoid_volume")
mesh.make("ellipsoid.msh")
```

| mesh | `r1` | `r2` | `mesh_size` | nodes | tetrahedrons | triangles | particles $N$ |
|------|------|------|-------------|-------|--------------|-----------|---------------|
| E1   | 4    | 8    | 0.5         | 3990  | 19339        | 3370      | 110795        |
| E2   | 20   | 30   | 2           | 5587  | 28096        | 4168      | 158571        |
| E3   | 60   | 80   | 3.5         | 22061 | 120211       | 11172     | 656632        |

E1 and E3 are reproduced exactly by this script with gmsh 4.15.2. The E2 file used here was
generated with an older gmsh: regenerated, it has 5662 nodes (1.3 % more), with the same adapted
height.

## Results

Laptop (Intel i7-8650U, 8 threads, WSL2), short runs ($10^{-11}$ s or $2 \times 10^{-12}$ s of
physical time), comparing `NbLevels = 6` (the code before commit 2fdda7690) with the dynamic
`NbLevels`. The magnetostatics per step is the time of
the two FMM calls of each step (potentials of $\mathbf{u}$ and $\mathbf{v}$, default
`FIRST_ORDER=OFF`), as printed in verbose mode (`magnetostatics done in ...`).

| mesh | dynamic `NbLevels` | magnetostatics per step, `NbLevels = 6` | magnetostatics per step, dynamic `NbLevels` | time per step, `NbLevels = 6` | time per step, dynamic `NbLevels` | speed-up |
|------|--------------------|-----------------------------------------|---------------------------------------------|-------------------------------|-----------------------------------|----------|
| E1   | 4      | 1130 ms                         | 169 ms                  | 1.34 s                | 0.35 s        | 3.8      |
| E2   | 4      | 2057 ms                         | 330 ms                  | 2.18 s                | 0.42 s        | 5.2      |
| E3   | 5      | 4214 ms                         | 1148 ms                 | 4.57 s                | 1.46 s        | 3.1      |

Accuracy: changing the height changes the FMM approximation, at the level of its own error. On
E1, the relative deviation of the total energy between `NbLevels = 6` and the dynamic `NbLevels`
after $10^{-11}$ s is
$2.5 \times 10^{-7}$; `ci-tests/full_test.py` (mesh `examples/ellipsoid.msh`, an ellipsoid of
semi-axes 4 and 8 nm with 167 nodes) gives an error of $10^{-6}$ on the magnetization and
$6.5 \times 10^{-7}$ on the energy (threshold $10^{-5}$).

## FMM algorithm built once

The FMM algorithm (`FmmClass`, i.e. `FFmmAlgorithmThreadTsm` of ScalFMM) is built once:
`fmm::algo` is a `std::unique_ptr<FmmClass>` created at the end of the constructor of
`scal_fmm::fmm`. Its constructor copies the kernel for each thread, with its precomputed rotation
tables, and reads the number of threads (`omp_get_max_threads()`): it is created after
`omp_set_num_threads(ScalfmmNbThreads)`, so that the FMM uses the number of threads of the
settings.

Before commit 29e2e6e35 ("Build the FMM algorithm once"), it was built at each call of
`fmm::demag`, i.e. twice per time step. Magnetostatics per step (two FMM calls), same laptop,
$2 \times 10^{-12}$ s of physical time, two runs each:

| mesh | algorithm built at each call | algorithm built once |
|------|------------------------------|----------------------|
| E1   | 178 to 202 ms                | 147 to 148 ms        |
| E2   | 408 to 418 ms                | 221 to 231 ms        |

The results are unchanged (relative deviations of $10^{-13}$, from the order of the parallel sums).

## Possible improvements

- The threshold was calibrated on one machine with 8 threads; it could be checked on a machine
  with many cores, where the near field (P2P) and the far field do not scale in the same way.
- A setting could override the computed height, for tests.
