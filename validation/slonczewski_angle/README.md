# Spin injection through a polarizer stack: angular validation on a thin free layer

`slonczewski_angle.py` checks the Robin boundary condition of the surface key `slonczewski`
(`src/slonczewskiInterface.h`, `spinAcc::stackBoundaryCondition`) in a non-collinear case, where
the mixing conductances act, against a reference written independently of feeLLGood.

The stack is F1 (ideal reservoir) / interface F1/N / spacer N of zero thickness / interface N/F2 /
F2 / reservoir ($s = 0$). F2 is the meshed thin layer (`thin_layer.geo`, 2 × 2 × 3 nm, h = 0.2 nm,
resolving $l_{sd}$ = 1 nm), with the magnetization $\mathbf{m} = \hat z$; the magnetization of F1
is $\mathbf{p} = (\sin\theta, 0, \cos\theta)$, $\theta$ = 0, 30, …, 180°. Both interfaces:
$AR^*$ = 0.5 fΩ·m², $\gamma$ = 0.75, $\mathrm{Re}\,G_{\uparrow\downarrow}$ = 0.55e15 S/m².

## Reference

- Interfaces: circuit theory of Brataas et al. (Phys. Rep. **427**, 157, 2006), in the variables
  of Valet–Fert ($\mu_s = (\mu_+-\mu_-)/2$): the spin current through an interface of
  magnetization $\mathbf{u}$ is
  ```math
  \mathbf{j}_s = \gamma\,j\,\mathbf{u} + \frac{1}{AR^*}\,\mathbf{u}\,(\mathbf{u}\cdot\Delta\boldsymbol{\mu}_s)/e
    + 2\,G_{\uparrow\downarrow}\,(\Delta\boldsymbol{\mu}_s)_\perp/e .
  ```
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
| 30  | 4.1525e-9 / 4.1504e-9 | 1.2808e-9 / 1.2940e-9 | 4.2e-3 |
| 60  | 6.5291e-9 / 6.5254e-9 | 2.0139e-9 / 2.0345e-9 | 1.1e-2 |
| 90  | 6.3616e-9 / 6.3578e-9 | 1.9622e-9 / 1.9823e-9 | 3.9e-3 |
| 120 | 4.4009e-9 / 4.3984e-9 | 1.3574e-9 / 1.3714e-9 | 1.7e-3 |
| 150 | 2.0581e-9 / 2.0569e-9 | 6.3480e-10 / 6.4132e-10 | 7.1e-4 |

At 0 and 180°, the transverse accumulation vanishes. The error decreases as $h^2$ (1.1e-2 at
h = 0.2 nm, 2.7e-3 at h = 0.1 nm). With the former conversion $k\,G_{\uparrow\downarrow}$, the error
reaches 170 %: the test fails.

The damping-like torque of this asymmetric structure follows the form of Slonczewski,
$\varepsilon(\theta)\sin\theta$ with $\varepsilon = C/[(\Lambda^2+1)+(\Lambda^2-1)\cos\theta]$, within
5 %, with $\Lambda^2 \approx 0.42$: larger torque near the parallel configuration. The ratio of the
field-like to the damping-like accumulation is independent of $\theta$ (about 0.31): it comes from
the precession of the transverse spin accumulation in F2, over $l_{sd}$.

Not covered: the sign of $\mathrm{Im}\,G_{\uparrow\downarrow}$ (zero here), and different densities
of states between the layers.
