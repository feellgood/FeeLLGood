# Tunnel barrier: 1D and angular validations

These scripts check the tunnel barrier of the surface key `tunnel` (`src/tunnelJunction.h`,
`electrostatSolver`, `spinAcc::tunnelBoundaryCondition`) with a uniform magnetization of the sample
F2, against references written independently of feeLLGood. The current concentration with a
non-uniform magnetization is checked by `validation/tunnel_two_domains`.

The barrier, between the electrode F1 (magnetization $\mathbf{p}$, potential $V_1$, not meshed) and
the sample F2 (magnetization $\mathbf{m}$ on the inlet), has the conductance per unit area of the
cosine law; the current leaving F2 through it and the spin current entering F2 are
```math
j_n = G_T(\mathbf{m})\,(V - V_1),\qquad
G_T = \frac{1 + c\,\mathbf{p}\cdot\mathbf{m}}{(1 + c)\,RA},\qquad
\mathbf{q} = \frac{\mu_B}{e}\,j_n\,\frac{P_1\,\mathbf{p} + P_2\,\mathbf{m}}{1 + c\,\mathbf{p}\cdot\mathbf{m}},
```
with $c = \mathrm{TMR}/(2 + \mathrm{TMR})$ and $P_1 P_2 = c$. Both cases use $RA$ = 1e-12 Ω·m²,
TMR = 1.5, $P_1$ = 0.7 ($P_2$ = 0.61), $V_1$ = −0.1 V (electrons flowing from F1 into F2), and
$\mathbf{m} = \hat z$, $\mathbf{p} = (\sin\theta, 0, \cos\theta)$.

## `wire_1d.py`: potential and longitudinal spin accumulation

Wire F2 (`wire.geo`, 4 × 4 × 100 nm, h = 1 nm), $\theta$ = 0, 90, 180°. The barrier is in series
with the wire, $j_n = -V_1/(1/G_T + L/\sigma)$, $V(0) = -j_n L/\sigma$. Along $\mathbf{m}$, the
problem is the one of `validation/valet_fert_1d` with the injected polarization
$P_l = (P_1\cos\theta + P_2)/(1 + c\cos\theta)$, solved analytically.

| θ (°) | $j_n$ (A/m²) | $V$(inlet), error | $s_z$(inlet) feeLLGood / ref (A/m) | max error on $s_z$ |
|-------|--------------|-------------------|------------------------------------|--------------------|
| 0   | 9.7561e10 | 4e-9  | 16.549 / 16.554 | 2.2e-3 |
| 90  | 6.8796e10 | 1e-8  | 5.891 / 5.893   | 2.2e-3 |
| 180 | 3.9604e10 | 1e-10 | −4.926 / −4.928 | 2.2e-3 |

The potential is exact (linear in $z$); the error on $s_z$ is the discretization error, as in
`validation/valet_fert_1d`. About 6 s.

## `thin_layer_angle.py`: transverse spin current and torques

Thin free layer F2 (`thin_layer.geo`, 2 × 2 × 3 nm, h = 0.2 nm, resolving $l_{sd}$ = 1 nm),
$\theta$ = 0, 30, …, 180°. The reference solves the equations of feeLLGood in F2 ($\beta = \beta'
= P$, anisotropic diffusion, relaxation, precession) exactly by eigenmodes, with the flux
$\mathbf{q}$ on the inlet and $\mathbf{s} = 0$ on the outlet, as in `validation/slonczewski_angle`.
The transverse part of $\mathbf{q}$, $\frac{\mu_B}{e}j_n P_1\sin\theta/(1 + c\cos\theta)$, is along
$x$; the torque on F2 is proportional to $\mathbf{m}\times\mathbf{s}$: $s_y$ gives the damping-like
torque, $s_x$ the field-like one.

| θ (°) | $j_n$ (A/m²) | damping-like $\int s_y\,dz$, feeLLGood / ref (A) | field-like $\int s_x\,dz$, feeLLGood / ref (A) | max error on $\mathbf{s}$ |
|-----|-----------|-------------------------|-------------------------|--------|
| 30  | 9.5912e10 | 6.7854e-10 / 6.7865e-10 | 1.4347e-10 / 1.4605e-10 | 1.6e-3 |
| 60  | 8.4946e10 | 1.1754e-9 / 1.1755e-9   | 2.4852e-10 / 2.5299e-10 | 4.1e-3 |
| 90  | 6.9963e10 | 1.3573e-9 / 1.3576e-9   | 2.8699e-10 / 2.9216e-10 | 1.2e-2 |
| 120 | 5.4977e10 | 1.1756e-9 / 1.1758e-9   | 2.4857e-10 / 2.5304e-10 | 1.1e-2 |
| 150 | 4.4005e10 | 6.7880e-10 / 6.7891e-10 | 1.4352e-10 / 1.4611e-10 | 2.7e-3 |

At 0 and 180°, the transverse accumulation vanishes. The error decreases as $h^2$ (1.2e-2 at
h = 0.2 nm, 2.9e-3 at h = 0.1 nm, at 90°; the field-like integral from 1.8 % to 0.5 %). The
tolerance is 2 % on the maximum error on $\mathbf{s}$ (about 14 s).

The torques being proportional to the transverse spin current entering F2, the damping-like torque
follows the form of Slonczewski, $\varepsilon(\theta)\sin\theta$ with
$\varepsilon = C/[(\Lambda^2+1)+(\Lambda^2-1)\cos\theta]$. Fitted to the integrals of feeLLGood:
- at fixed voltage ($\int s_y\,dz$): $\Lambda^2$ = 1.0004, the torque is nearly $\propto\sin\theta$,
  the barrier dominating the resistance;
- per unit current ($\int s_y\,dz/j_n$): $\Lambda^2$ = 2.5000 = 1 + TMR, larger torque near the
  antiparallel configuration (`dev-notes/slonczewski_polarisation.tex`, part "Jonction tunnel
  magnétique").

## Usage

```
./wire_1d.py [path to feellgood]            # default ../../src/feellgood
./thin_layer_angle.py [path to feellgood]
```
The meshes are built with `gmsh` if missing. The exit status is 0 if the errors are below the
tolerance at every angle.
