# Spin injection through a polarizer stack: angular validation on a thin free layer

`slonczewski_angle.py` checks the Robin boundary condition of the surface key `slonczewski`
(`src/slonczewskiInterface.h`, `spinAcc::stackBoundaryCondition`) in a non-collinear case, where
the mixing conductances act, against a reference written independently of feeLLGood.

The stack is F1 (ideal reservoir) / interface F1/N / spacer N of zero thickness / interface N/F2 /
F2 / reservoir ($s = 0$). F2 is the meshed thin layer (`thin_layer.geo`, 2 × 2 × 3 nm, h = 0.2 nm,
resolving $l_{sd}$ = 1 nm), with the magnetization $\mathbf{m} = \hat z$; the magnetization of F1
is $\mathbf{p} = (\sin\theta, 0, \cos\theta)$, $\theta$ = 0, 30, …, 180°. Both interfaces:
$AR^*$ = 0.5 fΩ·m², $\gamma$ = 0.75, $\mathrm{Re}\,G_{\uparrow\downarrow}$ = 0.55e15 S/m²;
$\mathrm{Im}\,G_{\uparrow\downarrow}$ = 0.15e15 (F1/N) and 0.25e15 S/m² (N/F2), larger than usual
so that its sign matters.

## Reference

- Interfaces: circuit theory of Brataas et al. (Phys. Rep. **427**, 157, 2006), in the variables
  of Valet–Fert ($\mu_s = (\mu_+-\mu_-)/2$, + being the majority spin, along the magnetization):
  the spin current through an interface of magnetization $\mathbf{u}$, from L to R, is
  ```math
  \mathbf{j}_s = \gamma\,j\,\mathbf{u} + \frac{1}{AR^*}\,\mathbf{u}\,(\mathbf{u}\cdot\Delta\boldsymbol{\mu}_s)/e
    + 2\,\mathrm{Re}\,G_{\uparrow\downarrow}\,(\Delta\boldsymbol{\mu}_s)_\perp/e
    + 2\,\mathrm{Im}\,G_{\uparrow\downarrow}\,\Delta\boldsymbol{\mu}_s\times\mathbf{u}/e ,
  \qquad \Delta\boldsymbol{\mu}_s = \boldsymbol{\mu}_{s,R} - \boldsymbol{\mu}_{s,L}.
  ```
  It is the equation of Brataas et al. (section on the circuit theory, the transverse spin current
  from N into F, $2\,\mathrm{Re}\,G_{\uparrow\downarrow}\,\mathbf{m}\times(\mathbf{V}_S\times\mathbf{m})
  + 2\,\mathrm{Im}\,G_{\uparrow\downarrow}\,\mathbf{V}_S\times\mathbf{m}$), with
  $G_{\uparrow\downarrow} = \frac{e^2}{h}\big[M - \sum r_\uparrow r_\downarrow^*\big]$ (Eq. 2,
  ↑ the majority spin, along $\mathbf{m}$, Fig. 2) and $\mathbf{V}_S = -\boldsymbol{\mu}_s/e$
  ($\mu = -eV$). The sign of the term in $\mathrm{Im}\,G_{\uparrow\downarrow}$ also follows from the
  scattering of a single channel: a transverse spin $(|{\uparrow}\rangle + |{\downarrow}\rangle)/\sqrt2$
  is reflected into $(r_\uparrow, r_\downarrow)/\sqrt2$, the transverse polarization absorbed is
  $1 - r_\uparrow^* r_\downarrow = \overline{G_{\uparrow\downarrow}}$ (per channel).
- F2: the equations of feeLLGood ($\beta = \beta' = P$, anisotropic diffusion, relaxation,
  precession), 1D, solved exactly by eigenmodes.
- $s = \mu_B N_0\,\mu_s/2$ (`validation/valet_fert_1d`).

This fixes the conversion of the mixing conductance: $G_{r} + iG_{i} = 2k\,G_{\uparrow\downarrow}$,
$k = 2/(e^2N_0)$, consistent with $G_\parallel = k/AR^*$: for a non-magnetic interface
($\gamma = 0$, $G_\uparrow = G_\downarrow = G_{\uparrow\downarrow}$), the conductance is isotropic.

## Usage

```
./slonczewski_angle.py [path to feellgood]    # default ../../src/feellgood
```
The mesh is built with `gmsh` if missing. The exit status is 0 if the maximum error on
$\mathbf{s}$, relative to $\max|\mathbf{s}|$, is below 2 % at every angle (about 15 s).

## Results (branch slonczewski)

The torque on F2 is proportional to $\mathbf{m}\times\mathbf{s}$: $s_y$ gives the damping-like
torque, $s_x$ the field-like one.

| θ (°) | damping-like $\int s_y\,dz$, feeLLGood / ref (A) | field-like $\int s_x\,dz$, feeLLGood / ref (A) | max error |
|-----|-----------------------|-----------------------|---------|
| 30  | 4.1054e-9 / 4.1033e-9 | 1.5626e-9 / 1.5755e-9 | 4.1e-3 |
| 60  | 6.6736e-9 / 6.6702e-9 | 2.4030e-9 / 2.4238e-9 | 1.1e-2 |
| 90  | 6.8434e-9 / 6.8403e-9 | 2.2136e-9 / 2.2346e-9 | 4.6e-3 |
| 120 | 5.0217e-9 / 5.0198e-9 | 1.3620e-9 / 1.3772e-9 | 1.9e-3 |
| 150 | 2.4700e-9 / 2.4694e-9 | 5.3402e-10 / 5.4143e-10 | 7.6e-4 |

At 0 and 180°, the transverse accumulation vanishes. The error decreases as $h^2$ (1.1e-2 at
h = 0.2 nm, 2.7e-3 at h = 0.1 nm, checked with $\mathrm{Im}\,G_{\uparrow\downarrow} = 0$). The
test fails if the conversion or the sign is wrong:
- with the former conversion $k\,G_{\uparrow\downarrow}$ instead of $2k\,G_{\uparrow\downarrow}$
  ($\mathrm{Im}\,G_{\uparrow\downarrow} = 0$): error 170 %;
- with the opposite sign of the term in $\mathrm{Im}\,G_{\uparrow\downarrow}$ in the reference
  (`IM_SIGN = -1`): error 13 %, the field-like accumulation differing by a factor 2 at 150°.

The damping-like torque of this asymmetric structure follows the form of Slonczewski,
$\varepsilon(\theta)\sin\theta$ with $\varepsilon = C/[(\Lambda^2+1)+(\Lambda^2-1)\cos\theta]$, within
5 %, with $\Lambda^2 \approx 0.54$ (0.42 with $\mathrm{Im}\,G_{\uparrow\downarrow} = 0$): larger
torque near the parallel configuration. Without $\mathrm{Im}\,G_{\uparrow\downarrow}$, the ratio of
the field-like to the damping-like accumulation is independent of $\theta$ (about 0.31): it comes
from the precession of the transverse spin accumulation in F2, over $l_{sd}$.

Not covered: different densities of states between the layers.
