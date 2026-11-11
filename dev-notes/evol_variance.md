# Variance of the magnetization in the `.evol` file

The `.evol` file records volume averages of the reduced magnetization,
`<Mx>`, `<My>` and `<Mz>`. Three more columns give the spatial spread of
each component around its average:

| column    | quantity                         |
|-----------|----------------------------------|
| `var(Mx)` | $\langle m_x^2\rangle - \langle m_x\rangle^2$ |
| `var(My)` | $\langle m_y^2\rangle - \langle m_y\rangle^2$ |
| `var(Mz)` | $\langle m_z^2\rangle - \langle m_z\rangle^2$ |

Like the averages, they can be restricted to a volume region with the
`region_name:` prefix:

```yaml
outputs:
  evol_columns: [t, <Mx>, <My>, <Mz>, var(Mx), var(My), var(Mz),
                 free_layer:<Mz>, free_layer:var(Mz)]
```

The variance is zero for a uniform state and grows when the magnetization
becomes non uniform (domain walls, vortices, spin waves). Its square root
is the standard deviation of the component over the volume.

## Definition

With $V$ the volume of the magnetic region (or of all the magnetic
regions) and $\mathbf{m} = \mathbf{M}/M_s$ the reduced magnetization:

```math
\langle m_i\rangle = \frac{1}{V}\int_V m_i \, dV, \qquad
\langle m_i^2\rangle = \frac{1}{V}\int_V m_i^2 \, dV, \qquad
\mathrm{var}(m_i) = \langle m_i^2\rangle - \langle m_i\rangle^2 .
```

Summing over the three components:

```math
\sum_i \mathrm{var}(m_i) = \langle |\mathbf{m}|^2\rangle - |\langle\mathbf{m}\rangle|^2 .
```

In the continuum $|\mathbf{m}| = 1$ and the sum is $1 - |\langle\mathbf{m}\rangle|^2$.
The discretized magnetization is only unitary at the nodes: inside a
tetrahedron the P1 interpolation of unit vectors has a norm smaller than
one, so that $\langle |\mathbf{m}|^2\rangle < 1$ and the sum is slightly
below $1 - |\langle\mathbf{m}\rangle|^2$. The gap is a measure of the
discretization error and vanishes when the mesh is refined.

## Variance restricted to a region

With the prefix, e.g. `free_layer:var(Mz)`, both averages are taken over
the region $R$ only: the tetrahedrons of other regions are skipped and the
integrals are divided by the volume $V_R$ of the region, not by the total
magnetic volume:

```math
\mathrm{var}_R(m_i) = \langle m_i^2\rangle_R - \langle m_i\rangle_R^2, \qquad
\langle f\rangle_R = \frac{1}{V_R}\int_R f \, dV .
```

It is the spread of $m_i$ inside the region, around the average of that
region. These variances are not additive. For a magnetic volume
$V = V_A + V_B$ split into two regions $A$ and $B$, the global variance is
the volume-weighted mean of the region variances plus the variance of the
region averages:

```math
\mathrm{var}(m_i) = \frac{V_A \, \mathrm{var}_A(m_i) + V_B \, \mathrm{var}_B(m_i)}{V}
+ \frac{V_A \left(\langle m_i\rangle_A - \langle m_i\rangle\right)^2
      + V_B \left(\langle m_i\rangle_B - \langle m_i\rangle\right)^2}{V} ,
```

with $\langle m_i\rangle = (V_A \langle m_i\rangle_A + V_B \langle m_i\rangle_B)/V$.
The formula extends to any number of regions. For example, a free layer
and a reference layer that are each uniformly magnetized, but in opposite
directions, have zero variance in each region while the global variance is
large: all of it comes from the second term.

## Discretization

On each tetrahedron $\mathbf{m}$ is interpolated linearly from its nodal
values $\mathbf{m}_k$ with the hat functions $a_k$:
$\mathbf{m} = \sum_k a_k \mathbf{m}_k$. Both integrals use the Gauss
quadrature of the element, as `<Mx>` did before:

```math
\int_{T} m_i \, dV \simeq \sum_{g} w_g \, m_i(\mathbf{x}_g), \qquad
\int_{T} m_i^2 \, dV \simeq \sum_{g} w_g \, m_i(\mathbf{x}_g)^2 .
```

With the default 5-point rule (exact up to degree 3), both are exact for
the P1 interpolation: $m_i$ is of degree 1 and $m_i^2$ of degree 2. The
reference values used by the unit test are

```math
\int_T a_k \, dV = \frac{V_T}{4}, \qquad
\int_T a_k a_l \, dV = \frac{V_T\,(1 + \delta_{kl})}{20},
```

hence

```math
\int_T m_i^2 \, dV = \frac{V_T}{20}\left(\sum_k m_{k,i}^2 + \Big(\sum_k m_{k,i}\Big)^2\right).
```

When feeLLGood is compiled with `ONE_GAUSS_POINT=ON`, $m_i^2$ is evaluated
at the barycenter only. This underestimates $\langle m_i^2\rangle$ (by
Jensen's inequality on each element), so the variance is biased low by a
quantity that goes to zero with the mesh size.

## Implementation

The computation is organized to read the nodal magnetization only once
per tetrahedron, whatever the number of requested columns.

1. `Tetra::Tet::integrals_u_u2()` (`src/tetra.h`) gathers the four nodal
   vectors in a $3\times 4$ matrix, interpolates them at the Gauss points
   with a single product by `eigen_a` ($3\times$`NPI`), and returns a
   $3\times 2$ matrix: first column $\int_T \mathbf{m} dV$, second column
   $\int_T (m_x^2, m_y^2, m_z^2) dV$. Both are two small matrix-vector
   products with the Gauss weights.
2. `Mesh::mesh::magMoments(region)` (`src/mesh.cpp`) sums these $3\times 2$
   matrices over the magnetic tetrahedrons with one parallel
   `std::transform_reduce` and divides by the volume. It returns a
   `MagMoments` structure holding the averages $\langle m_i\rangle$ and
   $\langle m_i^2\rangle$; `MagMoments::var()` forms the variances.
3. `Fem::saver()` (`src/save.cpp`) keeps a small cache, keyed by region
   index, valid for one output line. The first column that needs the
   moments of a region triggers the reduction; the other columns of the
   same region reuse the result. `<Mx>`, `<My>` and `<Mz>` use the same
   cache.

With the columns `<Mx> <My> <Mz> var(Mx) var(My) var(Mz)`, the mesh is
therefore traversed once instead of three times (previous `<Mi>` columns,
one `mesh::avg()` call each) or six times (naive implementation). The
cost of the reduction is dominated by the indirect reads of the nodal
magnetization; the extra squared terms are a negligible amount of
arithmetic. The reduction type is a fixed-size Eigen matrix, so the
reduction allocates nothing.

`mesh::avg()` is kept for `<dMx/dt>` and for the recentering, which work
on single components of $\mathbf{v}$ or $\mathbf{u}$.

## Numerical accuracy

The variance is computed as a difference of two averages. When it is
small compared with $\langle m_i\rangle^2$, the subtraction loses
relative precision: the absolute error is of the order of
$\varepsilon \langle m_i^2\rangle \le \varepsilon \approx 2\times 10^{-16}$.
This is far below any physically meaningful variance and below the
discretization error, so a two-pass (centered) or Welford-type algorithm
is not worth its cost. A rounding error can however make the difference
slightly negative for a uniform state; `MagMoments::var()` clamps the
result to zero, so that $\sqrt{\mathrm{var}}$ is always defined.

## Validation

- Unit test `ut_tetra/Tet_integrals_u_u2`: on a tetrahedron with random
  unit nodal vectors, `integrals_u_u2()` matches the exact P1 integrals
  above to rounding (absolute tolerance $10 \varepsilon V_T$).
- `ci-tests/full_test.py` still passes. The final $\langle\mathbf{m}\rangle$
  differs from the reference `ci-tests/full_test.evol` by $2\times 10^{-9}$,
  the same order as the difference between two runs of the current code
  ($3.5\times 10^{-9}$, parallel reductions).
- Same simulation with the new columns: the variances start at exactly
  zero (uniform initial state, clamped), then stay around $10^{-8}$; the
  columns with the `ellipsoid_volume:` prefix match the global ones to
  $6\times 10^{-17}$ (single region); the sum of the three variances is
  slightly below $1 - |\langle\mathbf{m}\rangle|^2$, as expected from the
  discussion above.
