# Spin injection through a polarizer stack: boundary conditions of the spin diffusion problem

This note explains the boundary conditions added to the spin diffusion problem of feeLLGood on the
branch `slonczewski`: where they come from, how and where they are implemented, what they can and
cannot describe, and a way to lift their main limitation; § 8 describes the injection through a
tunnel barrier of the branch `tunnel_junction`. It is written for readers with a
micromagnetic background but little experience of spintronics: the physical quantities are
introduced first. The derivations are detailed (in French) in
`dev-notes/slonczewski_polarisation.tex`.

## 1. Physical background

### 1.1 Spin accumulation, spin current, spin-transfer torque

In a metal, the conduction electrons carry a charge and a spin. In a ferromagnet, the electrons
whose spin is parallel to the local magnetization (majority electrons) and those whose spin is
antiparallel to it (minority electrons) conduct differently: a charge current flowing through a
ferromagnet is *spin polarized*, its spin is mostly along the magnetization.

When such a current enters another region (a non-magnetic metal, or a ferromagnet magnetized in
another direction), the spins it carries are not in equilibrium there: an excess of spin builds up
near the interface, the **spin accumulation**. It relaxes over two lengths:
- the **spin diffusion length** $`l_{sf}`$ (a few nm to hundreds of nm), over which spin-flip
  scattering brings the spins back to equilibrium;
- in a ferromagnet, the **transverse length** $`l_{sd}`$ (about 1 nm), over which the components of
  the spin accumulation perpendicular to the magnetization precess around the exchange field and
  dephase.

The angular momentum lost by the transverse part of the spin accumulation is transferred to the
magnetization: this is the **spin-transfer torque**, which can make a magnetization precess or
switch (spin valves, MRAM, spin-torque oscillators).

### 1.2 Quantities and parameters

| Symbol | Meaning | Unit | feeLLGood |
|---|---|---|---|
| $`\mathbf{m}`$ | unit vector of the magnetization (written $`\mathbf{u}`$ in the code) | – | solution of LLG |
| $`\mathbf{s}`$ | spin accumulation (a magnetization-like density, see § 4.4) | A/m | solution of the spin diffusion problem |
| $`\mathbf{J}`$, $`j_n`$ | charge current density, its normal component $`\mathbf{J}\cdot\mathbf{n}`$ on a surface, $`\mathbf{n}`$ being the **outward** normal ($`j_n > 0`$: the current leaves the sample, the electrons enter it) | A/m² | `jn` (number or `current_density`) |
| $`V`$ | electrostatic potential, $`\mathbf{J} = -\sigma\nabla V`$ | V | electrostatic solver |
| $`\sigma`$ | electrical conductivity | S/m | `sigma` |
| $`N_0`$ | density of states at the Fermi level | J⁻¹m⁻³ | `dens_state` |
| $`D_0`$ | spin diffusion constant, $`D_0 = 2\sigma/(e^2N_0)`$ | m²/s | computed |
| $`P`$ | bulk polarization rate of a ferromagnet, $`0 \le P < 1`$ ($`\beta = \beta' = P`$ of Zhang–Levy–Fert) | – | `P` |
| $`l_{sf}`$ | spin diffusion length (longitudinal in a ferromagnet) | m | `l_sf` |
| $`l_{sd}`$ | transverse (s–d exchange) length | m | `l_sd` |
| $`\tau_{sf}`$, $`\tau_{sd}`$ | spin-flip and precession times, $`1/\tau_{sf} = D_0(1-P^2)/l_{sf}^2`$ (ferromagnet), $`D_0/l_{sf}^2`$ (normal metal), $`1/\tau_{sd} = D_0/l_{sd}^2`$ | s | computed |
| $`\mu_B`$, $`e`$ | Bohr magneton, elementary charge ($`e > 0`$) | A·m², C | constants |
| $`\mathbf{p}`$ | magnetization of the reference layer F1 (unit vector) | – | `slonczewski.p` |
| $`AR^*`$, $`\gamma`$ | resistance–area product and spin asymmetry of an interface (Valet–Fert) | Ω·m², – | `AR`, `gamma` |
| $`G_{\uparrow\downarrow} = G_r + iG_i`$ | spin-mixing conductance of an interface (Brataas et al.), per unit area | S/m² | `Gr`, `Gi` |

Typical values (Co/Cu): $`AR^*`$ = 0.5e-15 Ω·m², $`\gamma`$ = 0.75,
$`\mathrm{Re}\,G_{\uparrow\downarrow}`$ about 0.5e15 S/m².

### 1.3 The spin diffusion problem solved by feeLLGood

The spin current is a tensor $`Q`$ of components $`Q_{i\alpha}`$: $`i`$ is the direction of flow
in space, $`\alpha`$ the component of the spin. Latin indices denote space, Greek indices the
spin, and the repeated indices are summed. In the model of Zhang, Levy and Fert, with
$`\beta = \beta' = P`$, it is
```math
Q_{i\alpha} = -P\,\frac{\mu_B}{e}\,J_i\,m_\alpha - D_0\left[\partial_i s_\alpha - P^2\,m_\alpha\,(\mathbf{m}\cdot\partial_i\mathbf{s})\right].
```
$`D_0 = 2\sigma/(e^2N_0)`$ is the diffusion constant of the code: it includes the factor 2 that
Zhang, Levy and Fert, and `dev-notes/slonczewski_polarisation.tex`, write in front of their
$`D_0`$ ($`2D_0`$ there); the note also writes the polarized term with
$`\mathbf{j}_e = \mathbf{J}`$ and $`|e|`$.
The first term is the spin carried by the charge current, polarized along $`\mathbf{m}`$ (it flows
with the electrons, opposite to $`\mathbf{J}`$); the second one is the diffusion of the spin
accumulation, with the diffusion constant $`D_0(1-P^2)`$ for the component of $`\mathbf{s}`$ along
$`\mathbf{m}`$ (longitudinal) and $`D_0`$ for the transverse ones. In a normal metal, $`P = 0`$.

The spin accumulation relaxes in picoseconds, while the magnetization evolves in nanoseconds:
feeLLGood solves the **steady state** at each time step, the magnetization being frozen,
```math
\partial_i Q_{i\alpha} + \frac{s_\alpha}{\tau_{sf}} + \frac{(\mathbf{s}\times\mathbf{m})_\alpha}{\tau_{sd}} = 0,
```
the last term being present in the ferromagnets only. The spin-transfer torque enters LLG as an
extra effective field proportional to the spin accumulation,
$`\mathbf{H}_{st} = \frac{D_0}{l_{sd}^2\,\gamma_0 M_s}\,\mathbf{s}`$ (`Tetra::calc_Hst`), i.e. a
torque $`\propto\mathbf{m}\times\mathbf{s}`$.

**Weak form.** With a test function $`\mathbf{w}`$ and $`\mathbf{q} = -\mathbf{n}\cdot Q`$
($`q_\alpha = -n_iQ_{i\alpha}`$) the spin current **entering** the sample through its boundary
$`\partial\Omega`$:
```math
\int_\Omega D_0\left[\partial_i s_\alpha - P^2 m_\alpha\,(\mathbf{m}\cdot\partial_i\mathbf{s})\right]\partial_i w_\alpha
+ \int_\Omega\left(\frac{\mathbf{s}}{\tau_{sf}} + \frac{\mathbf{s}\times\mathbf{m}}{\tau_{sd}}\right)\cdot\mathbf{w}
= \int_\Omega \frac{\mu_B}{e}P\sigma\,(\nabla V\cdot\nabla w_\alpha)\,m_\alpha
+ \int_{\partial\Omega}\mathbf{q}\cdot\mathbf{w}.
```
Code: `spinAcc::integrales` (diffusion, relaxation, precession, polarized current),
`Tet::add_spinDiffusionAniso` (the term in $`P^2`$), `spinAcc::solve` (assembly, boundary terms,
BiCGStab with the block preconditioner of `dev-notes/spin_accumulation_block_preconditioner.md`).
The unknowns are the three components of $`\mathbf{s}`$ at the nodes, with linear (P1) elements.

**Switching the anisotropic diffusion off.** The key `spin_accumulation.anisotropic_diffusion`
(default `true`) can be set to `false` to compare with the former model: the term in $`P^2`$ is then
dropped and the relaxation is $`1/\tau_{sf} = D_0/l_{sf}^2`$ in the ferromagnets too
($`\beta\beta' = 0`$, isotropic diffusion $`D_0`$). Both are switched together, so that $`l_{sf}`$
keeps the same meaning in both models (with the anisotropic term alone switched off, the
longitudinal decay length would become $`l_{sf}/\sqrt{1-P^2}`$). The polarizer stack of § 3 does not
depend on it; the equivalence with Valet–Fert (§ 4.4) and the validations of § 6 hold with `true`
only.

## 2. The boundary conditions before this branch

The boundary term $`\int_{\partial\Omega}\mathbf{q}\cdot\mathbf{w}`$ needs the spin current entering
through each surface. feeLLGood had three kinds of surfaces:
- **reservoir** (key `s`): $`\mathbf{s}`$ imposed (Dirichlet), usually 0, for a thick normal metal
  that absorbs the spin;
- **injection** (keys `jn` and `uP`): the entering spin current is imposed (Neumann),
  ```math
  \mathbf{q} = \frac{\mu_B}{e}\,j_n\,\mathbf{u}_P,
  ```
  $`\mathbf{u}_P`$ being a fixed unit vector describing a polarizer outside the mesh: `uP` is
  normalized when the settings are read, so that the injected spin current is always 100 %
  polarized;
- **free surface** (nothing given): $`\mathbf{q} = 0`$.

The injection condition is the limit of a very resistive polarizer. It ignores three effects that
govern the spin-transfer torque in a spin valve F1 / N / F2 (reference layer / non-magnetic spacer /
free layer):
1. **filtering by the interface of the free layer**: the spin carried by the charge current into
   F2 is polarized along the magnetization of F2 at the interface, $`\mathbf{m}_0`$, not along
   $`\mathbf{p}`$;
2. **backflow**: the spin accumulation built in F2 diffuses back into the spacer and changes the
   injected spin current;
3. **multiple reflections between F1 and F2**, which make the torque depend on the angle between
   $`\mathbf{m}_0`$ and $`\mathbf{p}`$ (Slonczewski's $`g(\theta)`$, or $`\varepsilon(\theta)`$ of
   mumax3 and OOMMF).

With $`\mathbf{q}`$ fixed, none of them can appear: $`\mathbf{q}`$ depends neither on
$`\mathbf{m}_0`$ nor on $`\mathbf{s}`$.

## 3. Origin of the new boundary condition

### 3.1 An interface in the circuit theory

The magnetoelectronic circuit theory (Brataas, Bauer, Kelly) describes an interface between a
normal metal N and a ferromagnet F of magnetization $`\mathbf{u}`$ by three conductances per unit
area: one for the spin along $`\mathbf{u}`$ (from $`AR^*`$ and $`\gamma`$) and the complex mixing
conductance $`G_{\uparrow\downarrow}`$ for the spin across $`\mathbf{u}`$. In the units of the spin
diffusion problem (§ 4.4), they form the 3×3 matrix
```math
\mathsf{G}(\mathbf{u}) = G_\parallel\,\mathbf{u}\mathbf{u}^T + G_r\,(\mathsf{I} - \mathbf{u}\mathbf{u}^T) - G_i\,[\mathbf{u}]_\times,
\qquad [\mathbf{u}]_\times\mathbf{w} = \mathbf{u}\times\mathbf{w}:
```
longitudinal conduction, absorption of the transverse spin (real part of the mixing conductance),
precession in the exchange field of the interface (imaginary part). The spin current from N into F
is then
```math
\mathbf{q} = P_I\,\frac{\mu_B}{e}\,j_n\,\mathbf{u} + \mathsf{G}(\mathbf{u})\,(\mathbf{s}_N - \mathbf{s}_F),
```
$`P_I = \gamma`$ being the polarization of the interface and $`\mathbf{s}_N`$, $`\mathbf{s}_F`$ the
spin accumulations on each side. This is linear in $`\mathbf{s}_F`$: a **Robin condition**, whose
coefficients depend on the magnetization.

### 3.2 Elimination of the spacer

The spacer N is thin (a few nm of Cu): its spin accumulation $`\mathbf{s}_N`$ is uniform across its
thickness and does not relax. The spin current entering N from F1 thus leaves it into F2, at each
point of the surface:
```math
P_1\frac{\mu_B}{e}j_n\,\mathbf{p} + \mathsf{G}_1(\mathbf{p})\,(\mathbf{s}_{F1} - \mathbf{s}_N)
= P_2\frac{\mu_B}{e}j_n\,\mathbf{m}_0 + \mathsf{G}_2(\mathbf{m}_0)\,(\mathbf{s}_N - \mathbf{s}_F),
```
indices 1 and 2 denoting the interfaces F1/N and N/F2. This 3×3 linear system gives
$`\mathbf{s}_N`$. Substituting it, with $`\mathsf{A} = \mathsf{G}_1 + \mathsf{G}_2`$, the spin
current entering F2 is
```math
\mathbf{q} = \frac{\mu_B}{e}\,j_n\,\mathbf{P}_{\mathrm{eff}}(\mathbf{m}_0) - \mathsf{G}_s\,\mathbf{s}_F,
\qquad
\mathbf{P}_{\mathrm{eff}} = \mathsf{G}_2\mathsf{A}^{-1}P_1\mathbf{p} + \mathsf{G}_1\mathsf{A}^{-1}P_2\mathbf{m}_0,
\qquad
\mathsf{G}_s = \mathsf{G}_2\mathsf{A}^{-1}\mathsf{G}_1,
```
with $`\mathbf{s}_{F1} = 0`$ (§ 3.3). The fixed vector $`\mathbf{u}_P`$ of the former condition is
replaced by an **effective polarization** $`\mathbf{P}_{\mathrm{eff}}`$ that interpolates between
the two polarizers, $`P_1\mathbf{p}`$ and $`P_2\mathbf{m}_0`$, with weights given by the
conductances, and the term $`-\mathsf{G}_s\mathbf{s}_F`$ is the backflow. The angular dependence
of the torque follows by itself, without any fitted function. Limits:
- F1/N very resistive ($`\mathsf{G}_1 \ll \mathsf{G}_2`$):
  $`\mathbf{P}_{\mathrm{eff}}\to P_1\mathbf{p}`$, $`\mathsf{G}_s\to 0`$, the former Neumann
  condition;
- collinear magnetizations: the series resistances of Valet–Fert;
- thin F2: the angular form $`\varepsilon(\theta)`$ of Slonczewski.

### 3.3 The bulk of F1

The interface F1/N needs the spin accumulation $`\mathbf{s}_{F1}`$ of F1 at the interface. F1
being uniformly magnetized along $`\mathbf{p}`$, its transverse spin accumulation vanishes within
about 1 nm, and its longitudinal one obeys a 1D collinear diffusion over its thickness $`t_1`$,
with a reservoir on its far side. Its solution is analytic: the bulk of F1 adds the spin
conductance
```math
g_F = \frac{2}{e^2N_0}\,\frac{\sigma_1(1-P_{F1}^2)}{l_{sf,1}}\coth\frac{t_1}{l_{sf,1}}
```
in series with the interface, which is equivalent to keeping $`\mathbf{s}_{F1} = 0`$ with the
renormalized parameters
```math
P_1^* = \frac{P_{F1}\,G_{\parallel1} + P_1\,g_F}{g_F + G_{\parallel1}},\qquad
G_{\parallel1}^* = \frac{g_F\,G_{\parallel1}}{g_F + G_{\parallel1}},
```
the mixing conductances being unchanged. Nothing has to be meshed in F1.

### 3.4 Transverse spin current through N/F2: `sink` or `series`

In the circuit theory, the ferromagnet absorbs the transverse spin current within a few atomic
layers: on the F side, the transverse spin accumulation is zero, and the transverse spin current
depends on $`\mathbf{s}_N`$ only. In the model of Zhang, Levy and Fert, the transverse spin
accumulation of F2 decays over $`l_{sd}`$ instead. Writing the mixing term with
$`\mathbf{s}_N - \mathbf{s}_F`$ would put the absorption by the interface and by the bulk of F2 in
series, and count it twice if $`G_{\uparrow\downarrow}`$ already describes it entirely (ab initio
values). The key `transverse` of `NF2` chooses:
- `sink` (default, circuit theory): $`\mathbf{s}_F`$ acts on N/F2 through its longitudinal part
  only, $`\mathsf{G}_s \to \mathsf{G}_s\,\mathbf{m}_0\mathbf{m}_0^T`$; the diffusion of F2 only
  spreads the absorbed spin current over $`l_{sd}`$;
- `series`: $`\mathbf{s}_N - \mathbf{s}_F`$; $`G_{\uparrow\downarrow}`$ should then be the part of
  the interface only (Petitjean et al.).

On the angular validation case, the damping-like torque differs by 15 to 17 % between them.

## 4. Implementation

### 4.1 Settings

The surface where the current enters the free layer replaces `uP` by a block `slonczewski`:
```yaml
inlet:
  jn: 1e12                       # or jn: current_density; jn > 0: electrons from F1 into F2
  slonczewski:
    p: [0, 0, 1]                 # magnetization of F1
    F1:  {sigma: 2e6, P: 0.5, l_sf: 20e-9, thickness: 15e-9}   # optional: bulk of F1
    F1N: {AR: 0.5e-15, gamma: 0.75, Gr: 5e14, Gi: 0}           # interface F1/N (required)
    NF2: {AR: 0.3e-15, gamma: 0.6, Gr: 5e14, Gi: 0,            # optional: interface N/F2
          transverse: sink}                                    # or series
```
- Without `F1`: $`\mathbf{s}_{F1} = 0`$ without renormalization (very conductive F1).
- Without `NF2` (spacer meshed, or perfect N/F2 interface):
  $`\mathbf{P}_{\mathrm{eff}} = P_1^*\mathbf{p}`$ and $`\mathsf{G}_s = \mathsf{G}_1(\mathbf{p})`$,
  independent of the magnetization.
- With `NF2`: the surface must border a magnetic region.
- `jn` is required on the same surface, `uP` and `s` are excluded. A time-dependent current density
  (`jn: current_density`, `dev-notes/current_density.md`) is accepted.

### 4.2 Where it is implemented

| What | Where |
|---|---|
| Reading and checking the block `slonczewski` | `src/settings.cpp`: `readStack`, `Settings::read` (surface regions), printed by `Settings::toYaml` |
| Stack in physical units, per surface | `Triangle::prm::slonczewski`, `Triangle::prm::stack` (`Slonczewski::StackSettings`), `src/triangle.h` |
| Conductance matrix, F1 renormalization, elimination of the spacer, unit conversion | `src/slonczewskiInterface.h`: `conductance`, `bulkConductance`, `renormalizeF1`, `effectiveBC`, `makeStack`, `Stack::bc` (pure functions, tested by `unit-tests/ut_slonczewskiInterface.cpp`) |
| Conversion with the density of states of the sample, checks | `spinAcc::makeStacks` (`src/spinAccumulationSolver.cpp`), once, at the construction of the solver |
| Boundary terms in the matrix and the right-hand side | `spinAcc::stackBoundaryCondition`, called by `spinAcc::solve` for every triangle of such a surface, at every solve |
| Anisotropic diffusion and longitudinal relaxation (needed for consistency with Valet–Fert), switched off with `anisotropic_diffusion: false` | `Tet::add_spinDiffusionAniso` (`src/tetra.cpp`), `spinAcc::integrales`, `spinAcc::anisotropicDiffusion`, `Settings::anisotropic_diffusion` |
| Warning about regions of different densities of states (§ 5.2) | `spinAcc::checkDensityOfStates` (verbose mode) |

### 4.3 Assembly

On each triangle of the surface, at each Gauss point:
1. the magnetization of F2, $`\mathbf{m}_0`$, is interpolated from the nodes (the latest one, as in
   the volume terms) and normalized;
2. `Stack::bc(m0)` computes $`\mathsf{G}_1(\mathbf{p})`$, $`\mathsf{G}_2(\mathbf{m}_0)`$, solves the
   3×3 system of the spacer (LU), and returns $`\mathbf{P}_{\mathrm{eff}}`$ and $`\mathsf{G}_s`$
   (times $`\mathbf{m}_0\mathbf{m}_0^T`$ with `sink`);
3. the boundary term $`\int_\Gamma\mathbf{q}\cdot\mathbf{w}`$ is split: the source
   $`\int_\Gamma\frac{\mu_B}{e}j_n\,\mathbf{P}_{\mathrm{eff}}\cdot\mathbf{w}`$ goes to the
   right-hand side, the term $`\int_\Gamma(\mathsf{G}_s\mathbf{s})\cdot\mathbf{w}`$ to the matrix (a
   boundary mass matrix with 3×3 blocks).

The matrix is assembled again at each solve, the magnetization having changed. The symmetric part
of $`\mathsf{G}_s`$ is positive and strengthens the coercivity; with `sink`, the boundary block is
not symmetric, which BiCGStab accepts. With a time-dependent current density, $`j_n`$ is multiplied
by `spinAcc::jn_factor` at the time of the solve.

### 4.4 Units: from the interface parameters to the spin diffusion problem

The interface parameters of the literature relate spin currents to the spin chemical potential
$`\mu_s = (\mu_+ - \mu_-)/2`$ of Valet–Fert ($`\mu_\pm`$: electrochemical potentials of the majority
and minority electrons). feeLLGood works with $`\mathbf{s}`$, in A/m:
```math
\mathbf{s} = \mu_B\,N_0\,\frac{\boldsymbol{\mu}_s}{2}.
```
With it, the bulk equations of Valet–Fert are exactly those of feeLLGood. The factor 1/2 comes from
$`D_0 = 2\sigma/(e^2N_0)`$, twice the coefficient of the Einstein relation: $`\mathbf{s}`$ is half
the physical magnetization $`\mu_BN_0\mu_s`$, which changes no result. An interface then gives, with
$`k = 2/(e^2N_0)`$,
```math
G_\parallel = \frac{k}{AR^*},\qquad P_I = \gamma,\qquad G_r + iG_i = 2k\,G_{\uparrow\downarrow}\quad[\mathrm{m/s}].
```
The factor 2 of the mixing part is checked on a non-magnetic interface ($`\gamma = 0`$), whose
conductance must be isotropic ($`G_\parallel = G_r`$, unit test `makeStack_isotropic`), and by the
angular validation, which fails with the factor 1. The stack not being meshed, all its conductances
(both interfaces and the bulk of F1) are converted with the single $`N_0`$ of the volume region
bordering the surface: `spinAcc::makeStacks` finds the tetrahedron of each triangle and stops if a
surface borders regions of different `dens_state`.

## 5. Limitations

### 5.1 Limitations of the model

- **One reference layer per surface, uniform and fixed**: $`\mathbf{p}`$ is a single vector, F1 is
  a 1D collinear layer folded into the interface. A non-uniform F1, a synthetic antiferromagnet or a
  free reference layer cannot be described.
- **Thin spacer**: uniform spin accumulation across it, no relaxation, a point-to-point relation
  between F1 and F2 (no lateral spreading in the spacer).
- **The stack is on the boundary of the mesh**: only the free layer F2 (and what lies beyond it) is
  meshed. The condition applies to an external surface, not to an interface inside the mesh.
- **Transverse length**: the transverse spin accumulation of F2 decays over $`l_{sd}`$ of about
  1 nm. With `sink`, the transverse spin current entering F2 is set on the surface whatever the
  mesh, but where it is absorbed over the thickness, hence the field-like torque, needs a mesh of
  the order of $`l_{sd}`$. On a free layer of 4 nm with the parameters of
  `examples/pillar_slonczewski_current_density_ramp1ns`, against the exact 1D solution (reference
  of `validation/slonczewski_angle`, θ = 90°, `sink`), a 2 nm mesh gives the damping-like torque
  within 3 % and overestimates the field-like torque by 29 % (see the README of the example).
- **Bulk torque only**: the whole spin current $`\mathbf{q}`$ enters F2 and the torque comes from
  $`\mathbf{m}\times\mathbf{s}`$ in the volume; a torque localized on the interface would need the
  spin currents on both sides of it (Petitjean et al.).
- **Charge transport**: the interface resistances $`AR^*`$ act on the spin only; the electrostatic
  potential is solved without them (no voltage drop at the interfaces, no magnetoresistance).

### 5.2 What the continuity of the spin accumulation imposes

The unknowns are the values of $`\mathbf{s}`$ at the nodes, with linear (P1) elements: a node shared
by two volume regions carries **one** value of $`\mathbf{s}`$, so that $`\mathbf{s}`$ is continuous
everywhere in the mesh, also across the boundary between two regions. The weak form also makes the
normal spin current continuous there (no interface term is assembled). Two consequences:

1. **No interface resistance between meshed regions.** At a real interface (F1/N, N/F2, Co/Cu…), the
   spin chemical potential **jumps**, the jump being set by the conductances $`AR^*`$, $`\gamma`$,
   $`G_{\uparrow\downarrow}`$ of § 3.1, and the interface itself polarizes the current ($`\gamma`$)
   and absorbs the transverse spin ($`G_r`$). A continuous $`\mathbf{s}`$ cannot represent any of
   this: meshing F1, N and F2 as three regions describes three bulk layers in perfect contact,
   without their interfaces, which usually dominate the spin transport in metallic multilayers. This
   is why the polarizer stack is not meshed but folded into a boundary condition: the jumps happen
   outside the mesh, on its boundary. The same limitation holds for the electrostatic potential.
2. **Different densities of states.** With $`\mathbf{s} = \mu_BN_0\boldsymbol{\mu}_s/2`$, the
   continuity of $`\mathbf{s}`$ between two regions A and B is
   ```math
   N_{0,A}\,\boldsymbol{\mu}_{s,A} = N_{0,B}\,\boldsymbol{\mu}_{s,B},
   ```
   not the continuity of $`\boldsymbol{\mu}_s`$ of Valet–Fert, unless the two regions have the same
   `dens_state`. On `validation/two_regions_1d` (a meshed Cu spacer and a ferromagnet with
   $`N_{0,F}/N_{0,N}`$ = 2.24), the spin accumulation of the ferromagnet at the interface differs by
   a factor about 2 between the two rules, and the torque with it. Workaround in the current code:
   give the same `dens_state` to all the regions, each one keeping its own `sigma`, hence its own
   $`D_0`$. In verbose mode, `spinAcc::checkDensityOfStates` warns about neighbouring regions of
   different `dens_state`.

## 6. Validation

| Case | What it checks | Result |
|---|---|---|
| uniform cylinder, injection along $`\mathbf{m}`$ | anisotropic diffusion, longitudinal $`l_{sf}`$ | decay length 10.07 nm for $`l_{sf}`$ = 10 nm, amplitude ratio $`1/(1+P)`$ exact |
| $`P_1 = 1`$, conductances nearly zero | limit of the former `uP` condition | same $`\mathbf{s}`$ within 5e-8 |
| `validation/valet_fert_1d` | collinear stack F1 / F1N / N / NF2 / F2 against the analytic Valet–Fert solution: conversion $`k/AR^*`$, renormalization of F1, signs | 0.26 % with $`h`$ = 1 nm, error in $`h^2`$; fails (8.5 %) with a wrong factor |
| `validation/slonczewski_angle` | thin F2, $`\theta`$ = 0 to 180°, mixing conductances and the sign of $`G_i`$, `sink` and `series`, against an independent reference | 1.1 % at $`h`$ = 0.2 nm, error in $`h^2`$; fails with a wrong conversion (170 %) or sign (13 %) |
| `validation/two_regions_1d` | continuity of $`\mathbf{s}`$ between regions of different `dens_state` | the model of feeLLGood within 0.14 % |
| `examples/pillar_slonczewski_current_density_ramp1ns` | spin-transfer switching of a perpendicular pillar under a ramp of current | switching AP → P, onset at the macrospin critical current density 3.5e11 A/m² |

## 7. A way around the continuity of the spin accumulation: discontinuous interfaces

The main limitation, § 5.2, comes from the single value of $`\mathbf{s}`$ at a node shared by two
regions. It can be lifted by letting $`\mathbf{s}`$ jump across selected internal interfaces, the
jump being governed by the interface law of the circuit theory, as is done for the thermal contact
resistance in heat transfer:

1. **Duplicate the nodes of the interface.** Each node of an internal surface $`\Gamma`$ between
   regions A and B (a physical surface of the mesh, declared as an interface in the settings) is
   split into a node of A and a node of B: $`\mathbf{s}`$ then has two traces, $`\mathbf{s}_A`$ and
   $`\mathbf{s}_B`$, on $`\Gamma`$. gmsh can produce such meshes (plugin `Crack`), or feeLLGood can
   duplicate the nodes when it reads the mesh.
2. **Assemble an interface term.** The spin current through $`\Gamma`$, from A into B, is the law of
   § 3.1 written between the two traces,
   ```math
   \mathbf{q}_\Gamma = P_I\,\frac{\mu_B}{e}\,j_n\,\mathbf{u} + \mathsf{G}(\mathbf{u})\,(\mathbf{s}_A - \mathbf{s}_B),
   ```
   and the weak form gains $`\int_\Gamma\mathbf{q}_\Gamma\cdot(\mathbf{w}_B - \mathbf{w}_A)`$: a
   matrix term coupling the two traces (the same 3×3 blocks as `stackBoundaryCondition`, with the
   signs of $`\mathbf{w}_A`$ and $`\mathbf{w}_B`$) and a source term. With the conductances
   converted for each side, the unknown can be $`\boldsymbol{\mu}_s`$ (or $`\mathbf{s}/N_0`$)
   instead of $`\mathbf{s}`$, which also gives the continuity of $`\boldsymbol{\mu}_s`$ of
   Valet–Fert between regions of different densities of states where no interface is declared.
3. **Do the same for the electrostatic potential**, with the interface resistance $`AR^*`$, so that
   the current distribution and the voltage drop include the interfaces.

The whole stack F1 / N / F2 could then be meshed (option 3 of the Slonczewski note, approach of
Abert et al.): non-uniform or free reference layers, synthetic antiferromagnets, several free
layers, lateral spreading in the spacer, and the interfacial torque from the difference of the spin
currents on both sides. The cost: the mesh must contain the thin spacer and the reference layer,
the nodes of each interface are doubled, and the boundary condition of § 3 remains useful as a
cheaper model when the stack is uniform. The present Robin condition is the special case where
everything on one side of $`\Gamma`$ is eliminated analytically.

## 8. Injection through a tunnel barrier (branch `tunnel_junction`)

The branch `tunnel_junction`, created from `slonczewski`, adds a second kind of injection: through
a tunnel barrier F1 / I (typically CoFeB/MgO), not meshed, key `tunnel` of a surface region. The
derivations and the discussion of the models are in the part "Jonction tunnel magnétique" of
`dev-notes/slonczewski_polarisation.tex`.

### 8.1 What changes compared with the spin valve

| | Spin valve (`slonczewski`) | Tunnel junction (`tunnel`) |
|---|---|---|
| spacer | diffusive metal, a node with a spin accumulation $`\mathbf{s}_N`$ | insulating barrier, coherent tunneling, no node |
| polarization | bulk and interfaces ($`\beta`$, $`\gamma`$, $`AR^*`$) | the barrier and its interfaces: effective tunnel polarizations |
| conductance | nearly independent of the angle (GMR of a few %) | cosine law between $`G_P`$ and $`G_{AP}`$, TMR of 100 to 600 % |
| feedback of $`\mathbf{s}`$ | $`AR^*`$ comparable to the spin resistance: Robin condition | $`RA`$ much larger than the spin resistance: flux condition |
| normal current density | imposed, uniform | depends on the local magnetization: Robin condition for the potential |

### 8.2 Conductance and spin current of the barrier

For any spin-conserving tunneling, the conductance per unit area of the barrier is linear in
$`\cos\theta = \mathbf{p}\cdot\mathbf{m}`$ (Slonczewski, 1989), $`\mathbf{p}`$ being the
magnetization of F1 and $`\mathbf{m}`$ the one of the sample on the surface:
```math
G_T(\mathbf{m}) = \frac{G_P + G_{AP}}{2} + \frac{G_P - G_{AP}}{2}\cos\theta
= \frac{1 + c\,\mathbf{p}\cdot\mathbf{m}}{(1 + c)\,RA_P},\qquad
c = \frac{G_P - G_{AP}}{G_P + G_{AP}} = \frac{\mathrm{TMR}}{2 + \mathrm{TMR}}.
```
This cosine law only needs the measured $`RA_P = 1/G_P`$ and TMR. The spin current is taken from
the Jullière model extended to non-collinear magnetizations: the spin current entering the sample
is $`(\mu_B/e)\,j_n\,\mathbf{P}_T`$ with
```math
\mathbf{P}_T(\mathbf{m}) = \frac{P_1\,\mathbf{p} + P_2\,\mathbf{m}}{1 + c\,\mathbf{p}\cdot\mathbf{m}},\qquad P_1P_2 = c,
```
which gives back the Jullière polarizations $`(P_1 + P_2)/(1 + P_1P_2)`$ (P) and
$`(P_2 - P_1)/(1 - P_1P_2)`$ (AP). The Jullière model is only used as a parametrization: $`P_1`$
and $`P_2`$ are effective tunnel polarizations, adjusted to the measured TMR, not the spin
polarizations of the densities of states of the electrodes (in Ni and Co, the latter is negative at
the Fermi level while the tunnel polarization is positive; with MgO, coherent tunneling through the
$`\Delta_1`$ states gives TMR far above the Jullière estimate). By default, $`P_1 = P_2 = \sqrt{c}`$
(symmetric junction): 0.65 for a TMR of 150 %.

The transverse part of the spin current, $`(\mu_B/e)\,j_n P_1\sin\theta/(1 + c\cos\theta)`$, gives
the damping-like torque. At fixed voltage, $`j_n \propto 1 + c\cos\theta`$ and the torque is
$`\propto\sin\theta`$ ($`\Lambda = 1`$ in the form of Slonczewski); per unit current,
$`\Lambda^2 = (1 + c)/(1 - c) = 1 + \mathrm{TMR}`$: the torque is larger near the antiparallel
state. These two results only rest on the cosine law.

### 8.3 Boundary conditions

- **Potential.** With $`V_1`$ the potential of F1 (equipotential), the current leaving the sample
  through the barrier is
  ```math
  j_n = G_T(\mathbf{m}_0)\,(V - V_1),
  ```
  a Robin condition: weak form $`\int_\Omega\sigma\nabla V\cdot\nabla w + \int_\Gamma G_T V w =
  \int_\Gamma G_T V_1 w`$. $`G_T`$ depending on the magnetization, the potential is solved again at
  each time step, and the current concentrates where $`\mathbf{m}`$ is parallel to $`\mathbf{p}`$.
- **Spin.** The feedback of the spin accumulation on the barrier is of the order of
  $`(k/RA)\,\mathbf{s}`$, against $`k/r_s`$ for the bulk of the sample,
  $`r_s = l_{sf}/(\sigma(1 - P^2))`$: with $`RA/r_s \sim 10^2`$ to $`10^4`$, it is neglected, and
  the condition is a flux,
  ```math
  \mathbf{q} = \frac{\mu_B}{e}\,j_n(\mathbf{x})\,\mathbf{P}_T(\mathbf{m}_0(\mathbf{x})).
  ```
  It is the former `uP` condition, with a polarization and a current density that both depend on
  the local magnetization.

### 8.4 Implementation

- `src/tunnelJunction.h`: pure functions `conductance` ($`G_T`$), `polarization` ($`\mathbf{P}_T`$),
  `polarizations` ($`P_1`$ by default, $`P_2 = c/P_1`$), `potentialF1` ($`V_1`$ or `voltage(t)`).
- `electrostatSolver`: on a `tunnel` surface, $`G_T(\mathbf{m}_0)\,a_ia_j`$ is added to the matrix
  and $`G_T(V_1 - V_D)\,a_i`$ to the right-hand side, $`\mathbf{m}_0`$ being interpolated and
  normalized at the Gauss points ($`V`$ is solved with $`V = 0`$ on the surface with `V`, $`V_D`$
  being added afterwards). The system is assembled again at each solve; the conjugate gradient starts
  from the previous $`V`$.
- `spinAcc`: the potential solver, formerly a local object of the constructor solved once, is a
  member; with a barrier, it is solved at the start of each `compute()`, before the spin.
  `tunnelBoundaryCondition` adds the flux $`\mathbf{q}`$ to the right-hand side, with
  $`j_n = G_T(\mathbf{m}_0)(V - V_1)`$ at the Gauss points.
- `compute_all` (`fem.h`) sets `voltage(t)` at each time step; column `voltage` of the `.evol` files.
- Unit tests `ut_tunnelJunction`: $`G_P = 1/RA`$, $`G_P/G_{AP} = 1 + \mathrm{TMR}`$, cosine law,
  default polarizations, Jullière limits, components of $`\mathbf{P}_T`$.

### 8.5 Settings

```yaml
voltage: "-0.8*t/1e-9"     # potential of F1 in V, function or expression of t (optional)
mesh:
  surface_regions:
    pillar_bottom:
      tunnel:
        p: [0, 0, -1]      # magnetization of F1
        RA: 0.5e-12        # Ohm m^2, parallel configuration
        TMR: 1.5           # (G_P - G_AP)/G_AP
        V: voltage         # or a number, in V
        P1: 0.65           # optional, sqrt(TMR/(2 + TMR)) by default
    pillar_top: {V: 0.0, s: [0, 0, 0]}
```
`tunnel` excludes `jn`, `V`, `s`, `uP` and `slonczewski` on the same surface, requires
$`c \le P_1 \le 1`$ (so that $`P_2 \le 1`$), and cannot be mixed with `jn: current_density`. With
$`V_1`$ below the potential of the sample, the electrons flow from F1 into the sample.

### 8.6 Validation and example

| Case | What it checks | Result |
|---|---|---|
| `validation/tunnel_junction/wire_1d.py` | wire of 100 nm, $`\theta`$ = 0, 90, 180°: potential and longitudinal spin accumulation against the analytic solution | potential exact (1e-8); $`s_z`$ within 2.2e-3 |
| `validation/tunnel_junction/thin_layer_angle.py` | thin layer of 3 nm, $`\theta`$ = 0 to 180°: transverse part and torques against eigenmodes of the sample | damping-like within 3e-4, field-like 1.8 %, max error on $`\mathbf{s}`$ 1.2e-2 in $`h^2`$; $`\Lambda^2`$ = 1.0004 at fixed voltage, 2.5000 = 1 + TMR per unit current |
| `validation/tunnel_two_domains` | slab with a Néel wall: current concentration and 2D potential against finite volumes | $`RA`$ = 1 Ω·µm²: $`j_P/j_{AP}`$ = 2.493, errors 6.1e-4 ($`V`$), 4.9e-5 ($`j_n`$); $`RA`$ = 2e-3 Ω·µm² (lateral spreading of 7.7 %): 5.0e-4, 8.7e-4 |
| `examples/pillar_tunnel_voltage_ramp1ns` | switching of a perpendicular MTJ pillar under a voltage ramp ($`RA_P`$ = 0.5 Ω·µm², TMR 150 %) | switching AP → P at −0.624 V (3.2 times the macrospin critical voltage 0.19 V); $`RA`$ from 1.248 to 0.505 Ω·µm²; current up to 1.9 times larger in the reversed part; solving the potential at each time step costs 3.7 % of a time step |

The example does not give the Oersted field (no `Bext`): it is small (about 2 mT at the edge of the
pillar at the critical voltage, 0.3 % of $`\mu_0H_{K,\mathrm{eff}}`$ = 0.73 T), azimuthal (zero
average over the disk), and the analytic field of a uniform current would be wrong, the current
being concentrated.

### 8.7 Open points

- Voltage dependence: TMR decreasing with $`V`$ and field-like torque in $`V^2`$ (Theodonis et al.,
  2006), not implemented.
- Current drive: the total current could be imposed using the linearity in $`V_1`$ at fixed
  magnetization.
- Oersted field computed from the current (Biot–Savart), a separate task.
- As for the spin valve, $`l_{sd}`$ is not resolved by a 2 nm mesh; the transverse spin current
  entering the sample is set by the boundary condition, only its distribution over the thickness
  depends on the mesh.

## 9. Correspondence with the notations of Baltz's textbook

The textbook of V. Baltz (chapters 3 and 4, exercise 6.3) introduces the physics this model rests
on: two-current model, Valet–Fert spin diffusion, interface resistances, spin transfer and circuit
theory. Its conventions differ from ours on three points:
- it takes $`e < 0`$;
- its arrows $`\uparrow, \downarrow`$ refer to the orientation of the *moment* ($`\uparrow`$ =
  majority, like the $`+`$ of Valet–Fert and Brataas), but it defines the spin accumulation and the
  spin current with respect to the *spin*, opposite to the moment:
  $`\mu_s^{\mathrm{B}} = \mu^\Uparrow - \mu^\Downarrow = -(\mu^\uparrow - \mu^\downarrow)`$ and
  $`J_s^{\mathrm{B}} = -(j^\uparrow - j^\downarrow)`$;
- its conductances are in S (not per unit area) and its torques in J/m³.

The superscript B denotes its quantities; $`k = 2/(e^2N_0)`$ as in § 4.4.

| Quantity | Baltz | feeLLGood | Relation |
|---|---|---|---|
| bulk polarization | $`\beta = (\alpha-1)/(\alpha+1)`$, $`\alpha = \rho^\downarrow/\rho^\uparrow`$ (§ 2.4) | `P` ($`\beta = \beta' = P`$) | $`P = \beta`$ |
| resistivities | $`\rho^{\uparrow(\downarrow)} = 2\rho^*(1\mp\beta)`$ | `sigma`, `P` | $`\rho^* = 1/\left(\sigma(1-P^2)\right)`$ |
| spin diffusion length | $`l_{sf}^*`$, $`1/l_{sf}^{*2} = 1/l_{sf}^{\uparrow2} + 1/l_{sf}^{\downarrow2}`$ (§ 3.1) | `l_sf` | $`l_{sf} = l_{sf}^*`$, $`1/\tau_{sf} = D_0(1-P^2)/l_{sf}^2`$ |
| spin resistance | $`\rho^* l_{sf}^*`$ in Ω·m² (§ 3.2) | $`g_F`$ (bulk of F1, § 3.3) | $`g_F = k\coth(t_1/l_{sf})/(\rho^*_{F1}\,l_{sf,1})`$ |
| interface | $`r^{\uparrow(\downarrow)} = 2r^*(1\mp\gamma)`$, eq. (40) | `AR`, `gamma` | $`AR^* = r^*`$, $`P_I = \gamma`$, $`G_\parallel = k/r^*`$ |
| spin memory loss | $`\delta`$ (§ 3.4) | – | not modelled |
| spin accumulation | $`\boldsymbol{\mu}_s^{\mathrm{B}}`$, along the spin | $`\mathbf{s}`$, $`\mu_s = (\mu_+ - \mu_-)/2`$ | $`s_\alpha = -\mu_BN_0\,\mu^{\mathrm{B}}_{s,\alpha}/4`$, i.e. $`\mu_s = -\mu_s^{\mathrm{B}}/2`$ |
| spin current | $`J_s^{\mathrm{B}}`$ in A/m², a tensor in 3D (§ 4.3) | $`Q_{i\alpha}`$ | $`Q_{i\alpha} = (\mu_B/\lvert e\rvert)\,J^{\mathrm{B}}_{s,i\alpha}`$ |
| angular momentum current | $`\mathcal{Q}^{\mathrm{B}} = \hbar J_s^{\mathrm{B}}/(2e)`$ | $`Q`$ | $`Q = -(2\mu_B/\hbar)\,\mathcal{Q}^{\mathrm{B}}`$ |
| transverse absorption | $`l_{sd} = \sqrt{\hbar D/J_{sd}}`$ (§ 4.3) | `l_sd`, $`1/\tau_{sd} = D_0/l_{sd}^2`$ | same role; compare the lengths, $`D`$ and $`J_{sd}`$ being defined differently |
| mixing conductance | $`G^{\uparrow\downarrow} = G_r^{\uparrow\downarrow} + iG_i^{\uparrow\downarrow}`$ in S; $`g^{\uparrow\downarrow} = AG^{\uparrow\downarrow}/(G_0/2)`$ in m⁻², $`G_0 = 2e^2/h`$ (§ 4.6) | `Gr`, `Gi` in S/m² | `Gr` + i `Gi` = $`G^{\uparrow\downarrow}/A = (e^2/h)\,g^{\uparrow\downarrow}`$; in the code $`G_r + iG_i = 2k\,(`$`Gr` + i `Gi`$`)`$ |
| torque | $`\mathbf{T}`$ in J/m³, $`\dot{\mathbf{M}} = -\lvert\gamma\rvert\mathbf{T}`$, eq. (65) | $`\boldsymbol{\tau} \propto \mathbf{m}\times\mathbf{s}`$ | – |

For example, $`g^{\uparrow\downarrow}`$ = 1e19 m⁻² gives `Gr` ≈ 3.9e14 S/m².

**Signs.** The relations for $`\mathbf{s}`$ and $`Q`$ are fixed by the physical results of the book
rather than by its intermediate equations. Far from the interface, in F,
$`J_s^{\mathrm{B}} = -\beta J_e`$ (eq. (34)), whereas here $`Q = -P(\mu_B/|e|)\,j_e\,\mathbf{m}`$:
hence $`Q = (\mu_B/|e|)\,J_s^{\mathrm{B}}`$. An excess of majority electrons gives
$`\mu_s^{\mathrm{B}} < 0`$ (p. 38) and $`\mathbf{s}\cdot\mathbf{m} > 0`$: hence the minus sign of
$`s_\alpha`$. Both relations agree with the diffusion term: with $`D_0 = 2\sigma/(e^2N_0)`$,
$`-D_0\nabla\mathbf{s}`$ gives back $`(\mu_B/|e|)\,J_s^{\mathrm{B}}`$ for Ohm's law
$`j = \sigma\nabla\mu/|e|`$. Eq. (22) of the book, $`j = \sigma\nabla\mu/e`$, would give the
opposite sign with $`e < 0`$: $`|e|`$ must be read there.

**Mixing conductance.** $`G^{\uparrow\downarrow}`$ has the same definition as in Brataas et al.
(eq. (80) of the book, $`\uparrow`$ = majority). The transverse current, eqs. (77)–(78), is linear
in the spin accumulation: going from the spin convention to the moment convention changes the sign
of both terms at once, without changing the ratio between the damping-like ($`G_r`$) and
field-like ($`G_i`$) torques. The sign of $`G_i`$ used here was checked against Brataas et al.
(`validation/slonczewski_angle`).

## References

- J. C. Slonczewski, J. Magn. Magn. Mater. **159**, L1 (1996); **247**, 324 (2002).
- T. Valet, A. Fert, Phys. Rev. B **48**, 7099 (1993): CPP giant magnetoresistance, collinear spin
  diffusion with interface resistances.
- S. Zhang, P. M. Levy, A. Fert, Phys. Rev. Lett. **88**, 236601 (2002): spin diffusion with a
  non-collinear magnetization.
- A. Brataas, G. E. W. Bauer, P. J. Kelly, Phys. Rep. **427**, 157 (2006): magnetoelectronic circuit
  theory, mixing conductance.
- C. Petitjean, D. Luc, X. Waintal, Phys. Rev. Lett. **109**, 117204 (2012): transverse spin
  currents in drift-diffusion, consistent interface conditions.
- C. Abert et al., Sci. Rep. **5**, 14855 (2015): 3D spin diffusion through the whole stack in
  micromagnetics.
- V. Baltz, *The Basics of Electron Transport in Spintronics*, EDP Sciences (2023): textbook,
  chapters 3 and 4; correspondence of the notations in § 9.
- M. Jullière, Phys. Lett. A **54**, 225 (1975); J. C. Slonczewski, Phys. Rev. B **39**, 6995
  (1989) and **71**, 024411 (2005): TMR, cosine law, torque in a tunnel junction.
- I. Theodonis et al., Phys. Rev. Lett. **97**, 237205 (2006): bias dependence of the torques.
- `dev-notes/slonczewski_polarisation.tex`: derivations, F1 options, conversion, validation, tunnel
  junction (in French).
