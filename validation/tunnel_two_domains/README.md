# Tunnel barrier: current concentration across a domain wall

`two_domains.py` checks the tunnel barrier of the surface key `tunnel` (`src/tunnelJunction.h`,
`electrostatSolver`, `spinAcc::tunnelBoundaryCondition`) with a non-uniform magnetization, against a
2D reference written independently of feeLLGood.

The sample F2 is a slab (`two_domains.geo`, 100 × 2 × 20 nm, h = 1 nm), the current flowing along
$z$: tunnel barrier on the inlet ($z = 0$), outlet at $V = 0$ ($z = t$). The magnetization of F1 is
$\mathbf{p} = \hat z$; the one of F2 is a Néel wall centred at $x = 0$,
$\mathbf{m} = (\sin\theta, 0, \cos\theta)$, $\theta = 2\arctan(e^{x/w})$, $w$ = 5 nm: parallel to
$\mathbf{p}$ on the left, antiparallel on the right. The conductance per unit area of the barrier
follows the cosine law and the current leaving F2 through it is
```math
j_n = G_T(\mathbf{m})\,(V - V_1),\qquad
G_T(\mathbf{m}) = \frac{1 + c\,\mathbf{p}\cdot\mathbf{m}}{(1 + c)\,RA},\qquad
c = \frac{\mathrm{TMR}}{2 + \mathrm{TMR}},
```
with TMR = 1.5 and $V_1$ = −0.1 V (electrons flowing from F1 into F2).

## Reference

The problem does not depend on $y$: $\nabla\cdot(\sigma\nabla V) = 0$ is solved in the plane
$(x, z)$ by second-order finite volumes on a grid of 0.2 nm (Neumann on $x = \pm L_x/2$, $V = 0$ on
$z = t$, the Robin condition above on $z = 0$), with the exact $\mathbf{m}(x)$. It needs `scipy`.
The grid error is about 1e-5 (0.2 nm against 0.1 nm).

Two cases:
- `local`: $RA$ = 1e-12 Ω·m², much larger than the resistance $t/\sigma$ = 5e-15 Ω·m² of F2: the
  current density is the one of a column, $j_n(x) = -V_1/(RA(\theta(x)) + t/\sigma)$, to 6e-4; far
  from the wall, $j_P/j_{AP}$ is close to 1 + TMR;
- `coupled`: $RA$ = 2e-15 Ω·m², comparable to $t/\sigma$: the current spreads laterally in F2 over
  about $t$ around the wall, 7.7 % from the local model of the columns. Unrealistic for a tunnel
  junction, this case checks the 2D solution of the potential, not only the Robin condition.

## Usage

```
./two_domains.py [path to feellgood]    # default ../../src/feellgood
```
The mesh is built with `gmsh` if missing. The exit status is 0 if the maximum errors on $V$ and
$j_n$ at the inlet, relative to their maximum, are below 1 % in both cases (about 10 s).

Note: the coordinates given to the function `initial_magnetization` are in m, whatever
`length_unit`.

## Results (branch tunnel_junction)

| case | mesh | error on $V$ | error on $j_n$ | $j_P/j_{AP}$ far from the wall: feeLLGood / columns |
|------|------|--------------|----------------|------------------------------------------------------|
| local   | h = 1 nm   | 6.1e-4 | 4.9e-5 | 2.4930 / 2.4925 |
| local   | h = 0.5 nm | 2.7e-4 | 1.2e-5 | |
| coupled | h = 1 nm   | 5.0e-4 | 8.7e-4 | 1.4412 / 1.4286 |
| coupled | h = 0.5 nm | 1.8e-4 | 2.9e-4 | |

In the local case, the current in the parallel domain is 2.49 times the one in the antiparallel
domain, close to 1 + TMR = 2.5, the resistance of F2 being small; the transition follows the wall.
In the coupled case, the ratio is smaller (series resistance of F2) and, far from the wall, still
1 % from the columns: the lateral spreading extends over a few $t$. The errors decrease with $h$,
as $h^2$ for $j_n$ in the local case.
