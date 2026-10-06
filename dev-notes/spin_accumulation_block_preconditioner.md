# Block preconditioner of the spin accumulation solver

The linear system of the spin accumulation is dominated, for short exchange lengths
$l_{sd}$, by the precession term $\frac{D_0}{l_{sd}^2}\,\mathbf{s}\times\mathbf{u}$, which couples
the three components of $\mathbf{s}$ at the same node. The diagonal (Jacobi) preconditioner ignores
this coupling, and the stabilised biconjugate gradient (BiCGStab) then needs hundreds of iterations.
Inverting exactly the $3\times3$ block of every node (block Jacobi preconditioner) removes it at a
negligible cost. This is the `BLOCK3` preconditioner of st-feeLLGood, implemented in feeLLGood by
`algebra::BlockDiagPrecond` (`src/algebra/block_precond.h`).

## The spin diffusion equation

The spin accumulation $\mathbf{s}$ (in A/m) is the solution of the steady-state equation
```math
-D_0\,\Delta\mathbf{s} + \frac{D_0}{l_{sf}^2}\,\mathbf{s}
  + \frac{D_0}{l_{sd}^2}\,\mathbf{s}\times\mathbf{u} = -\nabla\cdot J_s ,
\qquad D_0 = \frac{2\sigma}{e^2N_0},
\qquad (1)
```
where $\mathbf{u}$ is the reduced magnetization ($|\mathbf{u}| = 1$ in the ferromagnet,
$\mathbf{u} = 0$ elsewhere), $D_0$ the diffusion coefficient (conductivity $\sigma$, density of
states $N_0$), $l_{sf}$ the spin-flip length and $l_{sd}$ the s-d exchange length (`l_J` in
st-feeLLGood). The source $-\nabla\cdot J_s$ gathers the polarized current and the injection
through the faces; it only enters the right-hand side and plays no role here. The parameters are
constant per region.

## Discrete system

With P1 shape functions $a_i$ and test functions $\mathbf{w} = a_i\,\mathbf{e}_c$ ($c = x,y,z$),
the weak form of (1) reads (`spinAcc::integrales`, `src/spinAccumulationSolver.cpp`)
```math
\int D_0\,\nabla\mathbf{s}:\nabla\mathbf{w} + \int\frac{D_0}{l_{sf}^2}\,\mathbf{s}\cdot\mathbf{w}
 + \int\frac{D_0}{l_{sd}^2}\,(\mathbf{s}\times\mathbf{u})\cdot\mathbf{w} = \text{right-hand side},
```
the relaxation and precession terms being _lumped_ (mass concentrated on the nodes, with the nodal
value $\mathbf{u}_i$). Define, for every node $i$,
```math
S_{ij} = \int D_0\,\nabla a_i\cdot\nabla a_j,\qquad
\mu_i = \int\frac{D_0}{l_{sf}^2}\,a_i,\qquad
\nu_i = \int\frac{D_0}{l_{sd}^2}\,a_i ,
```
and the $3\times3$ matrix $C(\mathbf{u})$ of the map $\mathbf{s}\mapsto\mathbf{s}\times\mathbf{u}$:
```math
C(\mathbf{u}) = \begin{pmatrix} 0 & u_z & -u_y\\ -u_z & 0 & u_x\\ u_y & -u_x & 0\end{pmatrix},
\qquad C(\mathbf{u})\,\mathbf{s} = \mathbf{s}\times\mathbf{u} ,\qquad C(\mathbf{u})^T = -C(\mathbf{u}).
```
In feeLLGood, the unknowns are numbered node by node: component $c$ of node $i$ has the global
index $3i + c$ (st-feeLLGood numbers them by blocks of components, $cN + i$ for $N$ nodes). In
terms of nodal $3\times3$ blocks, the matrix of the system is
```math
K_{ij} = \big(S_{ij} + \delta_{ij}\,\mu_i\big)\,\mathrm{I}_3 + \delta_{ij}\,\nu_i\,C(\mathbf{u}_i),
\qquad (2)
```
that is, the sum of
- a _component-diagonal_ part (diffusion and relaxation), identical for the three components and
  coupling neighbouring nodes;
- a _node-diagonal, skew-symmetric_ part (precession), coupling the three components at the same
  node only.

$K$ is not symmetric, hence the use of BiCGStab.

## Why the diagonal preconditioner fails

The diagonal block of node $i$ is
```math
B_i = K_{ii} = d_i\,\mathrm{I}_3 + \nu_i\,C(\mathbf{u}_i),\qquad d_i = S_{ii} + \mu_i > 0 .
\qquad (3)
```
The diagonal of $K$ is $d_i$ (three times): the Jacobi preconditioner $D^{-1}$ keeps $d_i$ and
ignores $\nu_i\,C(\mathbf{u}_i)$. Since $C(\mathbf{u})\mathbf{u} = 0$ and
$C(\mathbf{u})^2 = \mathbf{u}\mathbf{u}^T - \mathrm{I}_3$ for $|\mathbf{u}| = 1$, the block $B_i$
has the eigenvalue $d_i$ (eigenvector $\mathbf{u}_i$) and the pair $d_i \pm \mathrm{i}\,\nu_i$
(plane orthogonal to $\mathbf{u}_i$). After the diagonal scaling, the local eigenvalues are $1$ and
$1 \pm \mathrm{i}\,\nu_i/d_i$.

The ratio $\nu_i/d_i$ compares the precession to the diffusion at the scale of the mesh. For a
mesh size $h$, $\int a_i \sim h^3$ and $S_{ii} \sim D_0\,h$, so that
```math
\frac{\nu_i}{d_i} \sim \frac{h^2/l_{sd}^2}{1 + h^2/l_{sf}^2} .
```
With $l_{sd} = 1$ nm and $h = 4$ nm, $\nu_i/d_i$ is of order 10: the diagonally scaled operator
has eigenvalues spread far along the imaginary axis, and BiCGStab converges slowly (about 500
iterations per solve, 1385 for the first one, on the mesh of the st-feeLLGood results below).

## Block Jacobi preconditioner

The preconditioner is the inverse of the block-diagonal part of $K$:
```math
M^{-1} = \operatorname{blockdiag}\big(B_1^{-1}, \dots, B_N^{-1}\big),
\qquad M^{-1} K = \mathrm{I} + \big(\text{coupling between nodes}\big).
```
It removes exactly the node-local skew-symmetric coupling; what remains is the diffusion between
neighbouring nodes, well conditioned at the scale of the mesh.

**Explicit inverse.** For $|\mathbf{u}_i| = 1$, the inverse of (3) is
```math
B_i^{-1} = \frac{d_i^2\,\mathrm{I}_3 - d_i\,\nu_i\,C(\mathbf{u}_i) + \nu_i^2\,\mathbf{u}_i\mathbf{u}_i^T}
                {d_i\,\big(d_i^2 + \nu_i^2\big)} ,
```
as checked with $C(\mathbf{u})\mathbf{u} = 0$ and
$C(\mathbf{u})^2 = \mathbf{u}\mathbf{u}^T - \mathrm{I}_3$:
$`(d\,\mathrm{I} + \nu C)(d^2\mathrm{I} - d\nu C + \nu^2\mathbf{u}\mathbf{u}^T) = d(d^2+\nu^2)\,\mathrm{I}`$.
In the non-magnetic regions ($\mathbf{u} = 0$), $B_i = d_i\,\mathrm{I}_3$ and the preconditioner
reduces to the diagonal one.

**Implementation.** The code does not use this formula: it reads the nine coefficients of each
block in the assembled matrix and inverts them (`Eigen::Matrix3d::inverse()`, by cofactors). It
thus stays exact whatever the content of the block (e.g. $|\mathbf{u}_i| \ne 1$, or additional
node-local terms in the future). The Dirichlet components (imposed spin accumulation `s` on some
faces) are handled as in the diagonal case: their rows and columns in the block are replaced by
those of the identity before the inversion, and the preconditioned vector is set to zero on them.
```
for each node i:
    a[r][c] = K(3i+r, 3i+c), r, c = 0..2   (identity for the Dirichlet components)
    inv_i   = a^{-1}
apply:  y[3i+c] = sum_c' inv_i[c][c'] * x[3i+c']
```
- `algebra::BlockDiagPrecond<BS>` (`src/algebra/block_precond.h`): construction (inverse of the
  `BS`x`BS` block of each node) and application `apply(x, y)`;
- `algebra::bicg_dir_prec` (`src/algebra/bicg.h`): BiCGStab with Dirichlet conditions and a general
  preconditioner; `algebra::bicg_dir` with Dirichlet values calls it with the diagonal
  preconditioner;
- `spinAcc::solve` (`src/spinAccumulationSolver.cpp`) builds `BlockDiagPrecond<3>` and calls
  `bicg_dir_prec`.

**Cost.** Nine coefficients per node are stored; the construction and each application are local
and cost a few tens of floating-point operations per node. Since $\mathbf{u}$ changes at every time
step, the blocks are rebuilt at every solve, at a cost negligible compared to one matrix-vector
product.

## Results

**st-feeLLGood** (measurements of the st-feeLLGood documentation). Cylinder of diameter 100 nm and
length 1000 nm, mesh size 4 nm (97 633 nodes, 534 095 tetrahedra), spin transport enabled,
6 accepted time steps, 8 threads. ILU(0) (incomplete LU factorization without fill-in,
recomputed at every solve) was also tried.

| preconditioner | diagonal | ILU(0) | $3\times3$ blocks |
|---|---|---|---|
| iterations, first solve | 1385 | 16 | 43 |
| iterations, following solves | about 500 | 7 | 14–16 |
| construction | – | 4 s | 20 ms |
| time of one solve (following) | about 70 s | about 6.3 s | about 2.2 s |
| spin accumulation (total) | 7 min 29 s | 59 s | 29 s |
| whole simulation | 8 min 56 s | 2 min 24 s | 1 min 54 s |

ILU(0) needs fewer iterations, but its factorization and its triangular solves are sequential and
cost more than the iterations saved. The block preconditioner captures precisely what makes the
system difficult (the precession, local to each node), at a negligible cost. The three
preconditioners give the same solution within the tolerance of the solver (relative deviation of
the outputs below 1e-5); on the small regression case `short_cyl_rich` of st-feeLLGood, the median
number of iterations drops from 310 (diagonal) to 43 (blocks) and 10 (ILU(0)).

**feeLLGood.**
- Domain wall in a cylinder of diameter 50 nm and length 1000 nm (15 003 nodes), $l_{sd} = 1$ nm,
  $l_{sf} = 10$ nm, $j_n = 10^{12}$ A/m², 28 time steps: 217 s with the diagonal preconditioner,
  43 s with the block one, same results (deviation of $\langle\mathbf{u}\rangle$ 4e-15).
- Unit test `spin_like_block_precond_solver` (`unit-tests/ut_algebra_bicg.cpp`), 2000 nodes with a
  strong skew-symmetric coupling of the components and Dirichlet conditions: about 500 to 800
  iterations with the diagonal preconditioner, 7 with the block one.

## Remarks

- The diffusion part of (2) is the same for the three components. An ILU or multigrid
  preconditioner of the scalar matrix $S + \operatorname{diag}(\mu)$, applied to each component,
  could be combined with the block one if the diffusion ever became the limiting factor (fine
  meshes with $h \ll l_{sd}$).
- When $h \ll l_{sd}$, $\nu_i/d_i \ll 1$ and the block and diagonal preconditioners become
  equivalent (the block one then brings no gain, at no noticeable cost).
