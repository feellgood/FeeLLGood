# pillar_slonczewski_current_density_ramp1ns: spin-transfer switching of a perpendicular pillar

Spin-transfer switching of the free layer of a perpendicular spin valve pillar, under a ramp of
the current density. The spin is injected through a polarizer stack that is not meshed (key
`slonczewski` of the bottom face): the spin current entering the free layer depends on its
magnetization (Robin boundary condition of the circuit theory), instead of being imposed with a
fixed polarization `uP`. It needs a feeLLGood built from the branch `slonczewski` (keys
`slonczewski` and `current_density`).

```bash
gmsh -3 pillar.geo -o pillar.msh       # built-in kernel of gmsh, as meshMaker.Cylinder
feellgood settings.yml
python3 stack_torque.py --evol pillar.evol
```

## The simulation

| | |
|---|---|
| Stack (bottom to top) | reservoir / F1, reference layer (10 nm, $\mathbf{p} = -\mathbf{e}_z$, not meshed) / interface F1/N / Cu spacer (not meshed) / interface N/F2 / **F2, free layer (meshed)** / Cu cap (reservoir) |
| Free layer F2 | disk of radius $R = 20$ nm and thickness $t = 4$ nm, axis $Oz$, $z \in [-2, 2]$ nm |
| Mesh | `pillar.msh`, size 2 nm: 1,117 nodes (`pillar.geo`, built-in kernel of gmsh) |
| Material of F2 | CoFeB-like with perpendicular anisotropy: $M_s = 1$ MA/m, $A = 2e-11$ J/m, $\alpha = 0.02$, $K = 8e5$ J/m³ along $z$ |
| Spin transport in F2 | $\sigma = 2e6$ S/m, `dens_state` 1.344e47 J⁻¹m⁻³, $P = 0.5$ ($\beta = \beta' = P$), $l_{sd} = 1$ nm, $l_{sf} = 5$ nm |
| Polarizer stack | F1: $\sigma = 2e6$ S/m, $P = 0.5$, $l_{sf} = 5$ nm, 10 nm; F1/N and N/F2: Co/Cu-like, $AR^* = 0.5$ fΩ·m², $\gamma = 0.75$, $G_{\uparrow\downarrow} = 0.55e15$ S/m² (Xia et al.), `transverse: sink` |
| Initial magnetization | close to $+\mathbf{e}_z$, antiparallel to $\mathbf{p}$, tilted by 5° towards $x$ |
| Current density | linear ramp `current_density: "2e12*t/1e-9"` on the bottom face: 0 at $t = 0$, 2e12 A/m² at 1 ns |
| Applied field | Oersted field of the current (`Bext`, space × time, calling `current_density(t)`) |
| Duration | 1 ns; one line of the `.evol` every ps, one `.sol` every 50 ps |

In `settings.yml`, the bottom face is:
```yaml
pillar_bottom:
  jn: current_density
  slonczewski:
    p: [0, 0, -1]
    F1: {sigma: 2.0e6, P: 0.5, l_sf: 5.0e-9, thickness: 10.0e-9}
    F1N: {AR: 0.5e-15, gamma: 0.75, Gr: 0.55e15, Gi: 0}
    NF2: {AR: 0.5e-15, gamma: 0.75, Gr: 0.55e15, Gi: 0, transverse: sink}
```
feeLLGood converts the stack once (verbose mode: `polarizer stack of pillar_bottom ... in m/s`)
with the `dens_state` of F2, $k = 2/(e^2N_0)$: $G_\parallel = k/AR^*$,
$G_r + iG_i = 2k\,G_{\uparrow\downarrow}$, the bulk of F1 being folded into F1/N
($P_1^* = 0.534$, $G_{\parallel1}^* = 1.56e5$ m/s), the spacer being eliminated. At each time
step, the spin current entering F2 through the bottom face is
```math
\mathbf{q} = \frac{\mu_B}{e}\,j_n\,\mathbf{P}_{\mathrm{eff}}(\mathbf{m}_0) - \mathsf{G}_s(\mathbf{m}_0)\,\mathbf{s},
```
$\mathbf{m}_0$ being the magnetization of F2 on the face, at its Gauss points. With
`transverse: sink`, the transverse spin current through N/F2 depends on the spin accumulation of
the spacer only (circuit theory): the mixing conductance describes the whole absorption, and the
diffusion in F2 spreads the absorbed spin current over $l_{sd}$. See the note
`slonczewski_polarisation.tex` and the validations `validation/valet_fert_1d`,
`validation/slonczewski_angle`.

**Sign.** `current_density(t)` is $\mathbf{j}\cdot\mathbf{n}$ on the bottom face, $\mathbf{n}$
being the outward normal $-\mathbf{e}_z$: `current_density(t)` > 0 is a current flowing along
$-z$, i.e. electrons flowing along $+z$, from F1 into F2. The spin-transfer torque then favours
$\mathbf{m}$ parallel to $\mathbf{p}$: F2, close to $+z$ (antiparallel to $\mathbf{p} = -z$),
switches towards $-z$ (AP → P) above the critical current density. The current density along $z$
is $j_z(t) = -$`current_density(t)`.

**Oersted field.** As in `cyl_oersted_current_density_ramp1ns`,
$\mathbf{B} = \frac{\mu_0 j_z}{2}(-y, x, 0)$ inside the pillar, given in `Bext` (feeLLGood does not
compute it): 0.025 T at the edge at the end of the ramp, small compared with the anisotropy field.

**Mesh and $l_{sd}$.** The transverse spin current is absorbed in F2 over $l_{sd} = 1$ nm, not
resolved by the 2 nm mesh. With `transverse: sink`, the transverse spin current entering F2 is set
by the circuit theory on the face, whatever the mesh; the mesh only acts on where it is absorbed
over the thickness, and on a small leak towards the reservoir of the top face. Compared with the
exact 1D solution (reference of `validation/slonczewski_angle`, F2 of 4 nm with the parameters of
this pillar, θ = 90°, `sink`), the integrals over the thickness of the transverse spin
accumulation are:

| mesh size | damping-like ($\int s_y$) | field-like ($\int s_x$) |
|---|---|---|
| 0.25 nm | −0.1 % | +0.6 % |
| 1 nm | −1.3 % | +11 % |
| 2 nm | −2.8 % | +29 % |
| 4 nm | −2.5 % | +40 % |

The damping-like torque, which drives the switching, is within 3 %; the field-like torque is
overestimated, but it is small here (about 7 % of the damping-like one). The micromagnetic
lengths are resolved: exchange length $\sqrt{2A/(\mu_0M_s^2)} = 5.6$ nm, wall width
$\sqrt{A/K_{\mathrm{eff}}} \approx 7.4$ nm.

**Computing time.** The mesh of 1 nm (6,063 nodes) makes LLG stiff (exchange field in $1/h^2$):
on 8 threads, 1.55 s per time step, failures of the LLG solver at `max(iter)`, time steps of about
1.5e-13 s, i.e. several hours for 1 ns. With 2 nm, the stiffness is divided by 4; with 4 nm (278
nodes), 0.12 s per time step, no failure, the time step reaching `max(dt)` = 5e-13 s. Whatever
the mesh, the spin accumulation costs 30 to 40 % of a time step, and the number of iterations of
LLG is the same with or without it.

## Macrospin estimate: `stack_torque.py`

The script converts the stack as feeLLGood does and gives, for a thin free layer whose spin
accumulation is neglected ($\mathbf{s}_F = 0$), the efficiency $\eta(\theta)$ of the damping-like
torque $(\hbar/2e)\,j_n\,\eta(\theta)\sin\theta$, $\theta$ being the angle between $\mathbf{m}$ and
$\mathbf{p}$:

| θ (°) | 1 | 30 | 60 | 90 | 120 | 150 | 179 |
|---|---|---|---|---|---|---|---|
| $\eta$ | 0.227 | 0.255 | 0.333 | 0.429 | 0.488 | 0.502 | 0.502 |

It has the asymmetric form of Slonczewski: the torque is larger near the antiparallel state.
The critical current density of the instability of the antiparallel state (linearized LLG,
perpendicular macrospin) is
```math
j_c = \frac{2e}{\hbar}\,\frac{\alpha\,\mu_0 M_s H_{K,\mathrm{eff}}\,t}{\eta(\pi)},\qquad
H_{K,\mathrm{eff}} = \frac{2K}{\mu_0 M_s} - (N_z - N_x)M_s,
```
with the demagnetizing factors of the disk $N_z = 0.797$, $N_x = N_y = 0.102$ (computed by the
script): $\mu_0 H_{K,\mathrm{eff}} = 0.727$ T, $\eta(\pi) = 0.502$, and
```math
j_c \approx 3.5e11\ \mathrm{A/m^2},
```
reached at 0.18 ns on the ramp. The field-like torque (zero here, $G_i = 0$ and
$\mathbf{s}_F = 0$), the Oersted field and the spin accumulation of F2 are neglected; the
simulation includes them, as well as the non-uniform magnetization.

With `--evol pillar.evol`, the script gives the current density at which $\langle M_z\rangle$
crosses 0 and plots the average magnetization against `current_density`
(`pillar_switching.png`).

## Results

Run of 2026-10-28 (branch `slonczewski`, commit `cdc8f86`), 8 threads: 22.5 min, 3,814 successful
time steps (mean 2.1e-13 s) and 4 failed ones, 0.35 s per time step.

The figure `pillar_switching.png` (`python3 stack_torque.py --evol pillar.evol`, not
committed) plots the average magnetization against `current_density`, the dashed line
being the macrospin critical current density.

| $t$ (ns) | `current_density` (A/m²) | $\langle M_x\rangle$ | $\langle M_y\rangle$ | $\langle M_z\rangle$ |
|---|---|---|---|---|
| 0.1 | 2.0e11 | +0.105 | −0.007 | +0.994 |
| 0.2 | 4.0e11 | +0.170 | −0.044 | +0.984 |
| 0.3 | 6.0e11 | +0.248 | −0.250 | +0.933 |
| 0.4 | 8.0e11 | −0.623 | −0.279 | +0.722 |
| 0.5 | 1.0e12 | −0.869 | −0.241 | +0.183 |
| 0.6 | 1.2e12 | −0.748 | +0.231 | −0.337 |
| 0.7 | 1.4e12 | −0.502 | −0.358 | −0.757 |
| 0.8 | 1.6e12 | +0.217 | −0.212 | −0.945 |
| 1.0 | 2.0e12 | +0.063 | −0.032 | −0.997 |

$\langle M_z\rangle$ goes below 0.9 at 0.328 ns (6.6e11 A/m²), crosses 0 at **0.536 ns,
1.07e12 A/m² (3.05 $j_c$)**, and goes below −0.9 at 0.760 ns (1.52e12 A/m²): F2 switches from
antiparallel to parallel to $\mathbf{p}$, as expected for electrons flowing from F1 into F2.

- **Below $j_c$**, the magnetization precesses around $+z$ (initial tilt of 5°) with an amplitude
  of about 0.1 that does not decay: the spin-transfer torque compensates part of the damping. The
  frequency, about 23 GHz (4 periods in 0.175 ns), agrees with
  $\gamma\mu_0 H_{K,\mathrm{eff}}/2\pi \approx 20$ GHz.
- **Above $j_c$**, the amplitude of the precession grows: this is the instability of the
  antiparallel state, which starts at the macrospin critical current density
  $j_c \approx 3.5e11$ A/m² of `stack_torque.py`.
- **Switching.** The ramp being fast (2e12 A/m² per ns, i.e. $j_c$ in 0.18 ns), the precession
  needs several periods to grow: $\langle M_z\rangle$ crosses 0 at 3 $j_c$. A slower ramp, or a
  constant current density, would give the switching time against the current density.
- **After switching**, the spin-transfer torque adds to the damping around $-z$: the precession is
  damped, and $\langle M_z\rangle = -0.997$ at 1 ns.

The time step goes down to about 5e-14 s during the switching (`max(du)` = 0.02). The output
gives a maximum angle of 41.3° between the magnetizations of neighbouring nodes: during the
switching, the magnetization is not uniform at the scale of the 2 nm mesh. A run with the 1 nm
mesh (several hours) would check the switching current density.
