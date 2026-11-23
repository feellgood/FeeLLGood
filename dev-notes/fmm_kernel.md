# Kernel of the FMM

## Context

At each time step, the magnetic scalar potential is computed at the magnetic nodes by the fast
multipole method (FMM) of ScalFMM: the potential of $\mathbf{u}$ (`phi`) and, with the default
`FIRST_ORDER=OFF`, the potential of $\mathbf{v}$ (`phiv`). The particles (magnetic nodes as
targets, Gauss points as sources) and the height of the octree are described in
`dev-notes/octree_height.md`.

The code is split in two parts:

- `src/fmm_demag.h` (class `scal_fmm::fmm`): what depends on the mesh, i.e. the positions of the
  particles, the height of the octree, the charges and the corrections of the triangles, for
  $\mathbf{u}$ and $\mathbf{v}$;
- `src/fmm_core.h` and `src/fmm_core.cpp` (class `scal_fmm::FmmCore`): the octree, the kernel and
  the algorithm of the FMM. The interface `fmm_core.h` depends neither on ScalFMM nor on Eigen:
  only `fmm_core.cpp` includes ScalFMM. Its headers are long to compile, and the BLAS and LAPACK
  declarations of `Utils/FBlas.hpp`, used by the Chebyshev kernel, conflict with those of the MKL
  headers, included by Eigen when feeLLGood is built with MKL (`EIGEN_USE_MKL_ALL`).

## Chebyshev kernel

The kernel is the symmetric Chebyshev kernel of ScalFMM, `FChebSymKernel`, with `ORDER = 6`
interpolation points per direction ($6^3 = 216$ per cell) and the interaction kernel $1/r$
(`FInterpMatrixKernelR`). Before commit 498cef3 ("Use the Chebyshev kernel for the FMM"), the
rotation kernel (spherical harmonics truncated at $P = 9$) was used.

- When the kernel is built, the M2L operators are compressed by singular value decompositions,
  using the symmetries of the interaction lists (a few ms); the M2L transfers are then matrix
  products. The kernel thus needs BLAS and LAPACK: CMake finds them, unless MKL is used, and
  ScalFMM must be built with BLAS support (`SCALFMM_USE_BLAS`, enabled when CMake finds a BLAS
  library while ScalFMM is configured). The CI installs OpenBLAS.
- The BLAS library runs on a single thread (`openblas_set_num_threads(1)` with OpenBLAS,
  `MKL_Set_Num_Threads(1)` with MKL, called through weak symbols): the FMM is already
  parallelized with OpenMP. With MKL, whose OpenMP runtime differs from the one of GCC, the BLAS
  threads would oversubscribe the cores so much that the run does not progress anymore.
- The kernel prints the time of its precomputations on `std::cout`: it is built with `std::cout`
  redirected, so that this line does not appear in the output of feeLLGood.

## Accuracy and computation time

The particles and the charges of $\mathbf{u}$ and $\mathbf{v}$ were dumped from feeLLGood during
simulations on the ellipsoids E1, E2, E3 of `dev-notes/octree_height.md` (after 100, 40 and 100
calls of the FMM), then the potentials were computed by a standalone program, with each kernel,
and compared to direct sums. Relative errors in norm 2 over the targets, time of the potentials
of $\mathbf{u}$ and $\mathbf{v}$ (two FMM passes, mean of 4 or 5 runs after a first one), Intel
i7-8650U, 8 threads, adapted height (4 for E1 and E2, 5 for E3):

| mesh | kernel        | time    | error (u) | error (v) |
|------|---------------|---------|-----------|-----------|
| E1   | rotation, P=9 | 148 ms  | 1.4e-5    | 4.4e-5    |
|      | Chebyshev, 4  | 109 ms  | 7.7e-5    | 2.0e-4    |
|      | Chebyshev, 5  | 114 ms  | 2.7e-5    | 1.5e-4    |
|      | Chebyshev, 6  | 122 ms  | 6.8e-7    | 3.1e-6    |
|      | Chebyshev, 7  | 146 ms  | 1.0e-7    | 5.8e-7    |
| E2   | rotation, P=9 | 210 ms  | 1.8e-5    | 2.7e-5    |
|      | Chebyshev, 4  | 147 ms  | 6.7e-5    | 1.1e-4    |
|      | Chebyshev, 5  | 148 ms  | 4.8e-5    | 7.8e-5    |
|      | Chebyshev, 6  | 163 ms  | 8.8e-7    | 1.8e-6    |
|      | Chebyshev, 7  | 203 ms  | 1.1e-7    | 2.2e-7    |
| E3   | rotation, P=9 | 919 ms  | 8.2e-4    | 5.2e-4    |
|      | Chebyshev, 4  | 495 ms  | 1.6e-3    | 1.1e-3    |
|      | Chebyshev, 5  | 465 ms  | 2.0e-3    | 1.6e-3    |
|      | Chebyshev, 6  | 636 ms  | 2.7e-5    | 2.2e-5    |
|      | Chebyshev, 7  | 1129 ms | 4.5e-6    | 3.4e-6    |

`ORDER = 6` is both faster and 20 to 30 times more accurate than the rotation kernel with $P = 9$.
`ORDER = 5` is barely faster and less accurate than the rotation kernel on E3; `ORDER = 7` is
slower, and needs larger OpenMP stacks (`OMP_STACKSIZE`) than the default ones.

The error of the rotation kernel on E3 does not depend on the height of the octree (8.0e-4,
8.2e-4 and 8.2e-4 for $h = 4, 5, 6$). The adapted height is also the fastest one with the
Chebyshev kernel (`ORDER = 6`, two passes):

| mesh | $h = 3$ | $h = 4$    | $h = 5$    | $h = 6$ |
|------|---------|------------|------------|---------|
| E1   | 444 ms  | **122 ms** | 154 ms     |         |
| E3   |         | 1950 ms    | **636 ms** | 3121 ms |

These measurements were made with ScalFMM1.5-lite, a lightened version of ScalFMM 1.5. With
ScalFMM 1.5 revision 22b9e4f6cf, the version installed by the CI, the times are similar (E1: 157 ms
with the rotation kernel, 130 ms with the Chebyshev kernel; E2: 211 ms and 183 ms), with the same
errors.

## Potentials of u and v in one pass

The potentials of $\mathbf{u}$ and $\mathbf{v}$ are sums over the same particles with different
charges: they are computed in a single FMM pass with two right-hand sides (`NVALS = 2`, or
`NVALS = 1` if `FIRST_ORDER` is ON). The cells, the particle containers and the kernel are
templated on `NVALS`; `FmmCore::compute` takes the charges and returns the potentials of each
right-hand side. Before commit 8730837 ("Compute the potentials of u and v in one FMM pass"), the
two potentials were computed by two passes.

The pass shares the traversals of the tree and, above all, the distances of the near field (P2P).
ScalFMM does not provide a suitable P2P for two right-hand sides: ScalFMM 1.5 has a scalar one,
which also computes the forces, and ScalFMM1.5-lite none. `src/fmm_p2p_two_rhs.h` specializes
`DirectInteractionComputer<double, 1, 2>` of ScalFMM: each distance is computed once for both
potentials, with the SIMD type of the P2P of ScalFMM (`SCALFMM_USE_AVX`, `SCALFMM_USE_AVX2`,
`SCALFMM_USE_SSE` or scalar). The target-source algorithm only calls its `P2PRemote`.

Same particles and machine as above, time of the potentials of $\mathbf{u}$ and $\mathbf{v}$, two
passes and one pass (the errors are those of the two passes):

| mesh | kernel        | two passes | one pass |
|------|---------------|------------|----------|
| E1   | Chebyshev, 5  | 114 ms     | 59 ms    |
|      | Chebyshev, 6  | 122 ms     | 72 ms    |
| E2   | Chebyshev, 5  | 148 ms     | 82 ms    |
|      | Chebyshev, 6  | 163 ms     | 101 ms   |
| E3   | Chebyshev, 5  | 465 ms     | 313 ms   |
|      | Chebyshev, 6  | 636 ms     | 490 ms   |

The adapted height is still the fastest one (`ORDER = 6`, one pass): E1 219, **72** and 140 ms for
$h = 3, 4, 5$; E3 1052, **490** and 3266 ms for $h = 4, 5, 6$.

## Results in feeLLGood

Magnetostatics per time step (the potentials of $\mathbf{u}$ and $\mathbf{v}$, as printed in
verbose mode), $2 \times 10^{-12}$ s of physical time, AMD EPYC 9654 (96 cores), comparing the
code before commit 47d3031 ("Isolate ScalFMM in its own translation unit"), i.e. the rotation
kernel with two passes, with the code after commit 8730837, i.e. the Chebyshev kernel
(`ORDER = 6`) with one pass. Both are compiled either with OpenBLAS or with MKL:

| mesh | threads | rotation, OpenBLAS | Chebyshev, OpenBLAS | rotation, MKL | Chebyshev, MKL |
|------|---------|--------------------|---------------------|---------------|----------------|
| E1   | 8       | 36.4 ms            | 19.0 ms             | 36.3 ms       | 21.0 ms        |
|      | 32      | 15.9 ms            | 9.9 ms              | 15.7 ms       | 10.5 ms        |
|      | 96      | 13.5 ms            | 8.3 ms              | 11.9 ms       | 9.4 ms         |
|      | 192     | 14.7 ms            | 9.8 ms              | 15.6 ms       | 10.7 ms        |
| E2   | 8       | 50.4 ms            | 27.9 ms             | 50.5 ms       | 31.5 ms        |
|      | 32      | 20.1 ms            | 14.7 ms             | 20.0 ms       | 15.4 ms        |
|      | 96      | 14.6 ms            | 12.2 ms             | 14.4 ms       | 13.2 ms        |
|      | 192     | 20.4 ms            | 13.1 ms             | 19.6 ms       | 13.7 ms        |
| E3   | 8       | 181.8 ms           | 107.8 ms            | 179.6 ms      | 139.3 ms       |
|      | 32      | 74.3 ms            | 52.5 ms             | 73.8 ms       | 64.3 ms        |
|      | 96      | 54.3 ms            | 48.0 ms             | 54.4 ms       | 50.1 ms        |
|      | 192     | 54.3 ms            | 45.3 ms             | 55.5 ms       | 49.4 ms        |

With 8 threads, the magnetostatics is 41 to 48 % faster with OpenBLAS, 22 to 42 % with MKL, whose
BLAS is slower than OpenBLAS for the small matrix products of the M2L transfers. With 32 to 192
threads, it is 8 to 39 % faster: the parts of `fmm::calc_demag` outside the FMM (charges,
corrections), which are sequential, weigh more. The total energy changes by
2e-9 to 4e-6 in relative value, at the level of the error of the rotation kernel, and
`ci-tests/full_test.py` passes with both libraries (deviation of 9.4e-7, threshold 1e-5).

## Possible improvements

- The charges and the corrections of the triangles are computed sequentially, for each right-hand
  side: parallelizing them would reduce the time with many threads.
- `ORDER = 7`, 5 to 8 times more accurate than `ORDER = 6`, would need larger OpenMP stacks.
