# Adiabatic and non-adiabatic torques on a wide domain wall

`wide_wall.py` checks the spin-transfer torque of a domain wall, computed by feeLLGood from the
spin accumulation and split by `tools/fg-torques` into its adiabatic and non-adiabatic parts,
against the local limit of Zhang–Li and its first non-local correction
(`dev-notes/adiabatic_torques.tex`). It checks the spin diffusion solver on a texture, with
$`\beta' = P'`$ = 0 or $`P`$ (key `P_prime`), and `fg-torques` itself.

The wire (`wire.geo`, 4 × 4 × 300 nm, h = 1 nm: P = 0.7, $`l_{sd}`$ = 1 nm, $`l_{sf}`$ = 5 nm)
carries a uniform current along $`z`$ (`jn` = 1e12 A/m² on the inlet, $`V = 0`$ and $`s = 0`$ on
the outlet). Its initial magnetization is a wall of width $`\delta`$ at its center $`z_0`$,
$`\theta(z) = 2\arctan e^{(z-z_0)/\delta}`$, $`\mathbf{m} = (\sin\theta, 0, \cos\theta)`$. The
magnetization is not relaxed: the spin accumulation of the iteration 0 is computed with this
profile, and `fg-torques` fits the form
$`\mathbf{T} = \eta\left[-(\mathbf{u}\cdot\nabla)\mathbf{m} + \beta\,\mathbf{m}\times(\mathbf{u}\cdot\nabla)\mathbf{m}\right]`$
to the torque $`-\mathbf{m}\times\mathbf{s}/(\tau_{sd}M_s)`$, by least squares weighted by
$`|(\mathbf{u}\cdot\nabla)\mathbf{m}|^2`$.

## Reference

In the local limit (texture wide compared with $`l_{sd}`$), with
$`\xi = \tau_{sd}/\tau_{sf} = (1-PP')\,l_{sd}^2/l_{sf}^2`$:
```math
\beta = \xi,\qquad \eta = \frac{1}{1+\xi^2}.
```
The first non-local correction comes from the transverse diffusion of the spin accumulation of the
local limit, $`l_{sd}^2\,\mathbf{s}_0''`$ in the source of the next order. For this wall, it changes
the local ratio of the two transverse components of $`\mathbf{s}`$ into
```math
\beta_{\mathrm{loc}}(z) = \xi - l_{sd}^2\,\frac{\theta'''}{\theta'}
= \xi - \frac{l_{sd}^2}{\delta^2}\left(1 - 2\,\mathrm{sech}^2\frac{z-z_0}{\delta}\right),
```
and the fitted $`\beta`$, average with the weight $`\theta'^2`$, into
```math
\beta_{\mathrm{fit}} = \xi + \frac{l_{sd}^2}{3\delta^2}.
```
The correction is absolute, the same for any $`P'`$; relative to $`\xi`$, it scales as
$`(l_{sf}/\delta)^2/(1-PP')`$: the wall must be wide compared with $`l_{sf}`$, not only $`l_{sd}`$,
for the local limit to hold.

## Usage

```
./wide_wall.py [path to feellgood]    # default ../../src/feellgood
```
The mesh is built with `gmsh` if missing (about 15 s). The exit status is 0 if, for
$`\delta`$ = 12.5 and 25 nm and $`P'`$ = 0 and 0.7, the fitted $`\beta`$ and the nodal
$`\beta_{\mathrm{loc}}`$ at the center of the wall match the first-order formulas within 5 % of
the non-local correction ($`l_{sd}^2/(3\delta^2)`$ and $`l_{sd}^2/\delta^2`$), and the fitted
$`\eta`$ matches $`1/(1+\xi^2)`$ within 2e-3.

## Results (branch tunnel_junction)

| $`\delta`$ | $`P'`$ | $`\xi`$ | correction / $`\xi`$ | $`\beta_{\mathrm{fit}}`$ feeLLGood / first order | error / correction | $`\beta_{\mathrm{loc}}(z_0)`$ feeLLGood / first order | error / correction | $`\eta_{\mathrm{fit}}`$ / $`1/(1+\xi^2)`$ |
|---|---|---|---|---|---|---|---|---|
| 12.5 nm | 0   | 0.0400 | 5.3 % | 0.042139 / 0.042133 | 0.3 % | 0.046455 / 0.046400 | 0.9 % | 0.99736 / 0.99840 |
| 12.5 nm | 0.7 | 0.0204 | 10 %  | 0.022568 / 0.022533 | 1.6 % | 0.026912 / 0.026800 | 1.7 % | 0.99877 / 0.99958 |
| 25 nm   | 0   | 0.0400 | 1.3 % | 0.040532 / 0.040533 | 0.2 % | 0.041621 / 0.041600 | 1.3 % | 0.99812 / 0.99840 |
| 25 nm   | 0.7 | 0.0204 | 2.6 % | 0.020936 / 0.020933 | 0.5 % | 0.022027 / 0.022000 | 1.7 % | 0.99935 / 0.99958 |

The largest errors on $`\beta`$ (δ = 12.5 nm, $`P'`$ = 0.7) are the next order, about
$`(l_{sd}/\delta)^4/\xi`$; the error on $`\eta`$ is a discretization error: with h = 2 nm, it is
three times larger (8.7e-4 instead of 2.9e-4 for δ = 25 nm, $`P'`$ = 0), and the error on the
correction of $`\beta`$ rises to 8 %.

The wire cannot be much longer with h = 1 nm: the electrostatic solver (conjugate gradient, 1000
iterations at most) does not converge on 600 layers. Hence $`\delta \le`$ 25 nm, $`z_0`$ being
6 $`\delta`$ at least from the ends.
