# pillar_tunnel_voltage_ramp1ns: spin-transfer switching of a perpendicular tunnel junction

Spin-transfer switching of the free layer of a perpendicular magnetic tunnel junction pillar, under
a ramp of the voltage. The current is injected through a tunnel barrier that is not meshed (key
`tunnel` of the bottom face): its conductance depends on the local magnetization of the free layer
(cosine law of the TMR), so that the current density is not uniform, and the potential is solved
again at each time step. It needs a feeLLGood built from the branch `tunnel_junction` (keys `tunnel`
and `voltage`). It is the tunnel counterpart of `pillar_slonczewski_current_density_ramp1ns`, with
the same free layer and mesh.

```bash
gmsh -3 pillar.geo -o pillar.msh       # built-in kernel of gmsh, as meshMaker.Cylinder
feellgood settings.yml
python3 junction_torque.py --evol pillar.evol
```

## The simulation

| | |
|---|---|
| Stack (bottom to top) | F1, reference layer ($\mathbf{p} = -\mathbf{e}_z$, equipotential at $V_1$, not meshed) / MgO barrier (not meshed) / **F2, free layer (meshed)** / cap (reservoir, $V = 0$) |
| Free layer F2 | disk of radius $R = 20$ nm and thickness $t = 4$ nm, axis $Oz$, $z \in [-2, 2]$ nm |
| Mesh | `pillar.msh`, size 2 nm: 1,117 nodes (`pillar.geo`, built-in kernel of gmsh) |
| Material of F2 | CoFeB-like with perpendicular anisotropy: $M_s = 1$ MA/m, $A = 2e-11$ J/m, $\alpha = 0.02$, $K = 8e5$ J/m³ along $z$ |
| Spin transport in F2 | $\sigma = 2e6$ S/m, `dens_state` 1.344e47 J⁻¹m⁻³, $P = 0.5$ ($\beta = \beta' = P$), $l_{sd} = 1$ nm, $l_{sf} = 5$ nm |
| Barrier | $RA_P = 0.5$ Ω·µm² (5e-13 Ω·m²), TMR = 150 % ($RA_{AP} = 1.25$ Ω·µm²), $P_1 = P_2 = 0.65$ (default, symmetric junction) |
| Initial magnetization | close to $+\mathbf{e}_z$, antiparallel to $\mathbf{p}$, tilted by 5° towards $x$ |
| Voltage | linear ramp `voltage: "-0.8*t/1e-9"`: potential of F1 from 0 at $t = 0$ to −0.8 V at 1 ns |
| Applied field | none (see *Oersted field* below) |
| Duration | 1 ns; one line of the `.evol` every ps, one `.sol` every 50 ps |

In `settings.yml`, the bottom face is:
```yaml
pillar_bottom:
  tunnel:
    p: [0, 0, -1]
    RA: 0.5e-12
    TMR: 1.5
    V: voltage
```
At each time step, the potential is solved with, on the bottom face, the Robin condition
```math
j_n = G_T(\mathbf{m}_0)\,(V - V_1),\qquad
G_T(\mathbf{m}_0) = \frac{1 + c\,\mathbf{p}\cdot\mathbf{m}_0}{(1 + c)\,RA_P},\qquad
c = \frac{\mathrm{TMR}}{2 + \mathrm{TMR}},
```
$\mathbf{m}_0$ being the magnetization of F2 on the face at its Gauss points, and the spin current
entering F2 is
```math
\mathbf{q} = \frac{\mu_B}{e}\,j_n\,\frac{P_1\,\mathbf{p} + P_2\,\mathbf{m}_0}{1 + c\,\mathbf{p}\cdot\mathbf{m}_0},
\qquad P_1P_2 = c,
```
the feedback of the spin accumulation of F2 on the barrier being neglected: $RA_P$ is 150 times
the spin resistance $l_{sf}/(\sigma(1 - P^2)) = 3.3\times10^{-15}$ Ω·m² of F2 (10³ to 10⁴ for
the usual $RA$ of 1 to 10 Ω·µm²). In F2, the spin current is the tensor $Q_{i\alpha}$ of the example
`pillar_slonczewski_current_density_ramp1ns`. See the part "Jonction tunnel magnétique" of
`dev-notes/slonczewski_polarisation.tex` and the validations `validation/tunnel_junction` and
`validation/tunnel_two_domains`.

**Sign.** With $V_1 < 0$, below the potential of F2, the current flows along $-z$ and the
electrons along $+z$, from F1 into F2. The spin-transfer torque then favours $\mathbf{m}$ parallel
to $\mathbf{p}$: F2, close to $+z$ (antiparallel to $\mathbf{p} = -z$), switches towards $-z$
(AP → P) above the critical voltage.

**Current concentration.** The conductance of the barrier is 2.5 times larger where $\mathbf{m}$ is
parallel to $\mathbf{p}$ than where it is antiparallel: during the switching, the current, hence
the torque, concentrates in the reversed parts of the free layer. At constant voltage, the current
increases by 2.5 when the junction switches from AP to P.

**Oersted field.** It is not given (no `Bext`), for three reasons:
- it is small: at the critical voltage, $j \approx 1.5\times10^{11}$ A/m², and
  $\mu_0 j R/2 \approx 2$ mT at the edge of the pillar, 0.3 % of the anisotropy field
  $\mu_0 H_{K,\mathrm{eff}} = 0.73$ T; at the end of the ramp, in the parallel state, it reaches
  20 mT, when the switching is over;
- it is azimuthal: its average over the disk vanishes, so that it does not act on the uniform
  magnetization, only on the curling of a non-uniform one;
- a `Bext` would be wrong: feeLLGood does not compute the field of the current, and the analytic
  field $\frac{\mu_0 j_z}{2}(-y, x, 0)$ of the spin valve example assumes a uniform current, which
  is not the case here (current concentration). Neglecting the field is more consistent than giving
  a wrong one.

**Mesh and $l_{sd}$.** As in `pillar_slonczewski_current_density_ramp1ns`, the 2 nm mesh does not
resolve $l_{sd} = 1$ nm; the transverse spin current entering F2 is set by the boundary condition,
the mesh only acts on where it is absorbed over the thickness (damping-like torque within 3 %, see
the README of that example).

**Computing time.** Solving the potential again at each time step costs 3.7 % of a time step: on
the first 20 ps, 0.386 s per time step with the barrier, 0.372 s with a fixed current density
($j_n$ and `uP`), for the same number of time steps (8 threads).

## Macrospin estimate: `junction_torque.py`

The transverse part of the spin current entering F2 gives the damping-like torque
$(\hbar/2e)\,j_n\,\eta(\theta)\sin\theta$, with the efficiency per unit current
$\eta(\theta) = P_1/(1 + c\cos\theta)$, $\theta$ being the angle between $\mathbf{m}$ and
$\mathbf{p}$:

| θ (°) | 0 | 30 | 60 | 90 | 120 | 150 | 180 |
|---|---|---|---|---|---|---|---|
| $\eta$ | 0.458 | 0.477 | 0.539 | 0.655 | 0.833 | 1.041 | 1.146 |

Per unit current, the torque is larger near the antiparallel state ($\Lambda^2 = 1 + \mathrm{TMR}$
in the form of Slonczewski); at fixed voltage, $j_n \propto 1 + c\cos\theta$ and the torque is
proportional to $\sin\theta$. The critical current density of the instability of the antiparallel
state (linearized LLG, perpendicular macrospin) and the critical voltage are
```math
j_c = \frac{2e}{\hbar}\,\frac{\alpha\,\mu_0 M_s H_{K,\mathrm{eff}}\,t}{\eta(\pi)},\qquad
|V_c| = j_c\,(RA_{AP} + t/\sigma),
```
$H_{K,\mathrm{eff}} = 2K/(\mu_0 M_s) - (N_z - N_x)M_s$, $N_z = 0.797$ and $N_x = 0.102$ being the
demagnetizing factors of the disk: $\mu_0 H_{K,\mathrm{eff}} = 0.73$ T,
$j_c = 1.54\times10^{11}$ A/m², $|V_c| = 0.19$ V, reached at 0.24 ns. Compared with the spin
valve example (efficiency 0.50 at θ = π), the efficiency of the barrier is 1.15: the critical
current density is 2.3 times smaller.

## Results

Branch `tunnel_junction`, 8 threads: 1 ns in 23 min (3,499 time steps, no failure, 0.40 s per
time step).

`junction_torque.py --evol pillar.evol` gives:
```
critical current density AP -> P (macrospin): jc = 1.542e+11 A/m^2, critical voltage |Vc| = 0.193 V
simulation: <Mz> crosses 0 at t = 7.8000e-10 s, voltage = -0.624 V (3.23 |Vc|)
```
F2 switches from AP to P: $\langle M_z\rangle$ goes from 0.996 to −0.982. Above $|V_c|$, the
magnetization precesses around $z$ with a growing amplitude, then reverses; $\langle M_z\rangle$
crosses 0 at 0.78 ns, at 3.2 times the macrospin critical voltage, as in the spin valve example
(3.1 times its critical current density): under a ramp, without thermal fluctuations, the
precession needs a voltage well above the threshold to grow within the ramp.

**Resistance.** The resistance-area product $-V_1/\langle j_n\rangle$, $\langle j_n\rangle$ being
the current density through the barrier averaged over the nodes of the bottom face (from the
`.sol` files, every 50 ps), goes from $RA_{AP}$ to $RA_P$ with the magnetization, as in a
measurement of the TMR:

| $t$ (ns) | 0.05 | 0.25 | 0.45 | 0.60 | 0.70 | 0.80 | 0.90 | 1.00 |
|---|---|---|---|---|---|---|---|---|
| $V_1$ (V) | −0.04 | −0.20 | −0.36 | −0.48 | −0.56 | −0.64 | −0.72 | −0.80 |
| $\langle j_n\rangle$ (A/m²) | 3.20e10 | 1.61e11 | 2.99e11 | 4.40e11 | 6.25e11 | 9.42e11 | 1.34e12 | 1.59e12 |
| $RA$ (Ω·µm²) | 1.248 | 1.243 | 1.206 | 1.091 | 0.896 | 0.680 | 0.536 | 0.505 |

($RA_{AP} + t/\sigma = 1.252$ and $RA_P + t/\sigma = 0.502$ Ω·µm²; the small differences come from
the average over the nodes and from the tilt of the magnetization.)

**Current concentration.** During the reversal, the magnetization of the bottom face is not
uniform, and the current density follows it (correlation of $j_n$ with $-m_z$: 1.000 at every
time):

| $t$ (ns) | $m_z$ on the bottom face | $j_n$ (A/m²) | max/min |
|---|---|---|---|
| 0.70 | 0.39 to 0.57 | 5.9e11 to 6.5e11 | 1.10 |
| 0.75 | −0.20 to 0.64 | 6.1e11 to 9.1e11 | 1.49 |
| 0.80 | −0.48 to 0.61 | 6.6e11 to 1.08e12 | 1.63 |
| 0.85 | −0.91 to 0.61 | 7.0e11 to 1.32e12 | 1.88 |
| 0.90 | −0.93 to −0.67 | 1.29e12 to 1.40e12 | 1.09 |

At 0.85 ns, the current density is 1.9 times larger in the reversed part than in the part not yet
reversed (2.5 times for fully parallel and antiparallel regions): the torque concentrates where
the reversal has started. This is the effect that a fixed current density (`jn`, `uP` or
`slonczewski`) cannot describe.

`junction_torque.py` also draws `pillar_switching.png`: the average magnetization and the
resistance-area product against the voltage.
