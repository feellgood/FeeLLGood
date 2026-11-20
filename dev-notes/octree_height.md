# Height of the FMM octree adapted to the particles

## Context

The demagnetizing field is computed in `src/fmm_demag.h` with the fast multipole method (FMM) of
ScalFMM 1.5 (rotation kernel, spherical harmonics truncated at $P = 9$). The particles are the
magnetic nodes (targets) and the Gauss points of the magnetic tetrahedrons and triangles (sources,
carrying the volume and surface magnetic charges). Their positions are scaled into the cube
$[-1, 1]^3$, inside a box of width 2.01 centered on the origin.

ScalFMM sorts the particles in an octree of height $h$: the leaves are the cells of depth
$h-1$, there are $2^{h-1}$ of them per direction. The potential at a target is the sum of

- a near field (P2P), computed directly from the particles of its leaf and of the neighbouring
  leaves;
- a far field, computed from multipole and local expansions (P2M, M2M, M2L, L2L, L2P) at all the
  levels of the octree.

Until now, the height was a constant `NbLevels = 6` whatever the mesh. On usual meshes (a few
thousand nodes, about $10^5$ particles) the leaves of depth 5 hold a few particles only: the far
field of these levels, and especially the M2L transfers (189 interactions per cell, each one
costing $O(P^3)$ with the rotation kernel), costs much more than the direct computation it
replaces.

## Choice of the height

The same problem was met and solved in st-feeLLGood (commit `6063be5`, `src/fmm_demag.cc`). The
heuristic is reproduced in `src/octree_height.h`:

> the height is the largest $h \in [3, 10]$ such that the non-empty leaves of the octree of
> height $h$ contain on average at least 150 particles.

The non-empty leaves are counted without building any octree:

1. the Morton key of each particle is computed at depth $9$ (the depth of the leaves of the
   deepest octree allowed, $2^9$ cells per direction): the bits of the three cell indices are
   interleaved, so that the key of the cell containing the particle at depth $h-1$ is the key
   shifted right by $3(10-h)$ bits;
2. the keys are sorted;
3. for $h = 10, 9, \dots, 4$, the number $L_h$ of distinct shifted keys is the number of
   non-empty leaves; the first $h$ such that $N \ge 150 L_h$ ($N$ number of particles) is
   returned, $h = 3$ otherwise.

Only the non-empty leaves are counted, so that the height adapts to the shape of the mesh: a thin
film or a long wire fills a small fraction of the cubic box, and gets a deeper octree than a
sphere with the same number of particles. The cost is a sort of the particles, done once when the
octree is built (a few ms for $10^5$ particles).

The number of sub-levels of the octree (`SizeSubLevels`, a parameter of the storage of the cells
in ScalFMM) is kept equal to 3, and reduced to $h-1$ if needed, as in st-feeLLGood.

### Calibration of the threshold

The threshold of 150 particles per leaf was calibrated in st-feeLLGood with the same kernel
(rotation, $P = 9$), on a short cylinder meshed with $2 \times 10^4$ to $7.5 \times 10^6$
particles: the fastest octree always has about 300 particles per non-empty leaf; one level deeper
(about 35 particles per leaf) is 4 to 5 times slower, one level less (about 2000 particles per
leaf) is 10 to 25 % slower. Since the number of particles per leaf is divided by about 8 from one
level to the next, the threshold 150 selects the fastest height.

| particles   | height | FMM time  | former height | former FMM time |
|-------------|--------|-----------|---------------|-----------------|
| 19614       | 3      | 7 ms      | 8             | 260 ms          |
| 1.3e5       | 4      | 96 ms     | 7             | 1.4 s           |
| 1e6         | 5      | 665 ms    | 7             | 7.9 s           |
| 7.5e6       | 6      | 6.1 s     | 7             | 26 s            |

(st-feeLLGood measurements, 8 threads.)

## Implementation in feeLLGood

- `src/octree_height.h` (no dependency on ScalFMM): `scal_fmm::mortonKey` and
  `scal_fmm::octreeHeight`, with the constants `MinLevels = 3`, `MaxLevels = 10` and
  `MinParticlesPerLeaf = 150`.
- `src/fmm_demag.h`: the constant `NbLevels` is removed. The class `scal_fmm::fmm` gets a member
  `nbLevels`, initialized from the normalized positions of the particles
  (`fmm::particlePositions`, same order and same scaling as the insertion in the tree) before the
  tree and the kernel, which are built with this height. The normalization coefficient `norm` is
  now initialized in the initializer list (it is needed by `particlePositions`).
- `src/main.cpp`: in verbose mode, the height is printed with the insertion of the particles:
  `Magnetostatics: particles inserted, using 8 threads, octree of 4 levels, in 28.21 ms`.
- `unit-tests/ut_octree_height.cpp`: Morton keys of the 8 cells of depth 1, uniform grid of
  $64^3$ particles filling the box (height 4), the same grid in one cell of depth 4 (height 8:
  the empty leaves are ignored), few particles (height 3), all particles at the same point
  (height 10).

## Results

Laptop (Intel i7-8650U, 8 threads, WSL2), short runs ($10^{-11}$ s or $2 \times 10^{-12}$ s of
physical time) of existing simulations, compared to master (`NbLevels = 6`). The time per step
includes the two FMM calls of each step (potentials of $\mathbf{u}$ and $\mathbf{v}$, default
`FIRST_ORDER=OFF`).

| mesh                                    | nodes | particles | height | magnetostatics per step, master | magnetostatics per step | time per step, master | time per step | speed-up |
|-----------------------------------------|-------|-----------|--------|---------------------------------|-------------------------|-----------------------|---------------|----------|
| ellipsoid, example0.1                   | 3990  | 110795    | 4      | 1130 ms                         | 169 ms                  | 1.34 s                | 0.35 s        | 3.8      |
| ellipsoid, fmr_5587nodes                | 5587  | 158571    | 4      | 2057 ms                         | 330 ms                  | 2.18 s                | 0.42 s        | 5.2      |
| ellipsoid, example0_enhanced_curling_v2 | 22061 | 656632    | 5      | 4214 ms                         | 1148 ms                 | 4.57 s                | 1.46 s        | 3.1      |

Accuracy: changing the height changes the FMM approximation, at the level of its own error. On
example0.1, the relative deviation of the total energy from master after $10^{-11}$ s is
$2.5 \times 10^{-7}$; `ci-tests/full_test.py` gives an error of $10^{-6}$ on the magnetization
and $6.5 \times 10^{-7}$ on the energy (threshold $10^{-5}$).

## FMM algorithm built once

The FMM algorithm (`FmmClass`, i.e. `FFmmAlgorithmThreadTsm` of ScalFMM) was built at each call of
`fmm::demag`, i.e. twice per time step. Its constructor copies the kernel for each thread, with its
precomputed rotation tables, and reads the number of threads (`omp_get_max_threads()`). It is now
built once, as in st-feeLLGood: `fmm::algo` is a `std::unique_ptr<FmmClass>` created at the end of
the constructor of `scal_fmm::fmm`, after `omp_set_num_threads(ScalfmmNbThreads)`, so that the
number of threads is unchanged.

Magnetostatics per step (two FMM calls), same machine, $2 \times 10^{-12}$ s of physical time, two
runs each:

| mesh                     | algorithm built at each call | algorithm built once |
|--------------------------|------------------------------|----------------------|
| ellipsoid, example0.1    | 178 to 202 ms                | 147 to 148 ms        |
| ellipsoid, fmr_5587nodes | 408 to 418 ms                | 221 to 231 ms        |

The results are unchanged (relative deviations of $10^{-13}$, from the order of the parallel sums).

## Possible improvements

- The threshold was calibrated on one machine and one kernel; it could be checked on a machine
  with many cores, where the near field (P2P) and the far field do not scale in the same way.
- A setting could override the computed height, for tests.
