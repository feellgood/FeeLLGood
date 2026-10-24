# Time-dependent current density

The key `current_density` of the settings gives the normal current density $j_n(t)$, in A/m^2, of
the surface regions that declare `jn: current_density`, $t$ being the time in seconds. It replaces
the constant value of `jn` of those surfaces. In st-feeLLGood, the same current is
`pulse(t)` × `J_n`, `pulse(t)` being a dimensionless factor defined in the case library
(`pulse.cc`). Here it is a JavaScript function given in the settings, like `Bext`, and the
functions of `Bext` can call it.

## Usage

### Definition

The surface region where the current is injected declares `jn: current_density`, and the function
is given at the top level of the settings, either as a JavaScript function of `t`:
```yaml
current_density: "function(t) { return t < 1e-9 ? 0 : -1e12; }"

mesh:
  surface_regions:
    inlet:
      jn: current_density
      uP: [0, 1, 0]
    outlet:
      V: 0
      s: [0, 0, 0]
```
or as an expression of `t`, which feeLLGood wraps into `function(t) { return (...); }`:
```yaml
current_density: "-1e12*min(t/1e-10, 1)"
```
The functions of the Math object are available without the `Math.` prefix (`min`, `exp`, `sin`,
`PI`...), as in `Bext`. A numeric `jn` is still accepted for a constant current density.

The shape is arbitrary (step, ramp, trapezoid, Gaussian, sinusoid...) and so is the sign.
Example of a trapezoid, the multi-line form of YAML keeping the function readable:
```yaml
current_density: |
  function(t) {
    const j0 = -1e12;                                                  // [A/m^2]
    const t_on = 1e-10, rise = 5e-11, length = 1e-9, fall = 5e-11;  // [s]
    if (t < t_on) return 0;
    if (t < t_on + rise) return j0*(t - t_on)/rise;
    if (t < t_on + rise + length) return j0;
    if (t < t_on + rise + length + fall) return j0*(1 - (t - t_on - rise - length)/fall);
    return 0;
  }
```
and of a sinusoidal current of 1 GHz switched on with a linear rise over 0.5 ns:
```yaml
current_density: "-1e12*min(t/5e-10, 1)*sin(2*PI*1e9*t)"
```

### Sign convention

As `jn`, `current_density(t)` is $\mathbf{j}\cdot\mathbf{n}$ with $\mathbf{n}$ the **outward**
normal of the surface: a positive value makes the current leave the sample through that surface.
To inject a current through the face $z = z_{\min}$ of a wire along $z$ (outward normal
$-\mathbf{e}_z$), so that it flows along $+z$, `current_density(t)` is negative:
$j_z(t) = -$`current_density(t)`.

### Effect on the current

With `spin_accumulation.enable: true`, the normal current density of the surface with
`jn: current_density` is `current_density(t)` at time $t$. Without spin accumulation the
electric problem is not solved, and `current_density` only acts through the applied field (see
below).

### Applied field following the current

The three forms of `Bext` can call `current_density(t)`, whether the spin accumulation is enabled
or not:
```yaml
# function of t
Bext: "function(t) { return [1e-13*current_density(t), 0, 0]; }"

# three expressions of t
Bext: ["1e-13*current_density(t)", 0, 0]

# product of a space function and a time function: only the time function can call current_density
Bext:
  space: "function(x, y, z) { return [-y, x, 0]; }"
  time: "function(t) { return -2*PI*1e-7*current_density(t); }"
```
The last form is the natural one for the Oersted field: its space profile is computed once on the
Gauss points, and its amplitude follows the current. In a wire along $z$ carrying the uniform
current density $j_z$, the Oersted field is $`\mathbf{B} = \frac{\mu_0 j_z}{2}(-y, x, 0)`$, in T
with $x$, $y$ in m; with the injection through the face $z = z_{\min}$,
$j_z = -$`current_density(t)` and $\mu_0/2 = 2\pi\times 1e-7$, hence the time function above.
feeLLGood does not compute the Oersted field of the current itself: the profile is the user's
responsibility.

If `Bext` calls `current_density` while the key `current_density` is not defined, the run stops
at the first evaluation of the field with the message
`Script ReferenceError: identifier 'current_density' undefined`.

### Checks

With the spin accumulation enabled:
- a surface region with `jn: current_density` while the key `current_density` is not defined stops
  the run: `Error: a surface region has jn: current_density, but the function current_density is
  not defined.`
- the key `current_density` defined while no surface region has `jn: current_density` gives a
  warning on std::cerr: the function then only acts through `Bext`.

The electrostatic problem requires at least one surface with `jn` and a single surface with `V`.
Several surfaces may have `jn`, either all numeric (each with its value) or all
`jn: current_density` (the same function for all of them); mixing them stops the run:
`Error: incorrect boundary conditions for potential V solver: the surfaces with jn cannot mix
jn: current_density and numeric jn.`

### Outputs

- `.evol`: column `current_density`, to be added to `outputs.evol_columns`, value of
  `current_density(t)` [A/m^2] at the time of the line (NaN if it is not defined). The columns
  `Hx`, `Hy`, `Hz` give the applied field including the current (form `Bext` function of `t` or
  three expressions only).
- `.sol` (with the spin accumulation): the columns are `idx mx my mz phi sx sy sz V`. The spin
  accumulation is the one of the time of the file, thus with the current of that time, and the
  last column `V` is the electrostatic potential at that time,
  $`V(t) = V_D + (j_n(t)/j_{\text{ref}})V_{\text{ref}}`$ (see the hypotheses), in V; with a
  numeric `jn` it is $V_0$. It costs no solve: $`V_{\text{ref}}`$ (or $V_0$) is kept from the
  start. When a `.sol` file is used as the
  initial magnetization (restart), the columns after `phi` are skipped: the spin accumulation and
  the potential are computed again at the start.
- `spin_accumulation.V_file`: the potential solved at the start, written once in
  `<file_basename>_V.sol` of the current directory: $V_0$ for a numeric `jn`, the reference
  potential $`V_{\text{ref}}`$ for `jn: current_density`, i.e. with jn = +1 A/m^2 on those surfaces
  and V = 0 on the surface with `V` (see the hypotheses). For a current injected through the
  surface ($j_n < 0$), $`V_{\text{ref}}`$ has the opposite sign of the actual potential.
- At start-up, the settings printed in verbose mode show the key `current_density`, and `J:
  current_density` for the surface region.

### From Python

`settingsMaker` has no dedicated method: set the keys of the dictionary directly,
```python
mySettings["current_density"] = "-1e12*min(t/1e-10, 1)"
mySettings["mesh"]["surface_regions"]["inlet"]["jn"] = "current_density"
```

## Hypotheses

1. **A single time profile.** One or several surfaces may have `jn: current_density`: the same
   function `current_density(t)` is then their normal current density, so that the potential is
   proportional to the single factor $`j_n(t)/j_{\text{ref}}`$. Surfaces with independent
   profiles, or a mix of numeric `jn` and `jn: current_density`, would need one reference solve
   per profile (superposition) and are not possible. The electric problem has a single surface
   with `V`.

2. **Linear and quasi-static electrokinetics.** The potential is the solution of
   ```math
   \nabla\cdot(\sigma\nabla V) = 0,\qquad
   -\sigma\frac{\partial V}{\partial n} = j_n(t) \text{ on the surface with } j_n,\qquad
   V = V_D \text{ on the surface with } V.
   ```
   The conductivity $\sigma$ does not depend on the magnetization (no anisotropic
   magnetoresistance), on the current or on the temperature (no Joule heating), and the current
   follows $j_n(t)$ instantaneously (no capacitive or inductive effect). The problem being linear,
   ```math
   V(t) = V_D + j_n(t)\,G,\qquad G = \frac{V_{\text{ref}}}{j_{\text{ref}}},
   ```
   where $`V_{\text{ref}}`$ is the solution of the reference problem: the reference current
   density $`j_{\text{ref}} = 1`$ A/m^2 on the surfaces with `jn: current_density` and **V = 0** on
   the surface with `V`. $G$, in V m^2/A, is the potential per unit normal current density.
   $`V_{\text{ref}}`$ is solved once, at the start, and never solved again. The factor
   $`j_n(t)/j_{\text{ref}}`$ is dimensionless, and numerically equal to the value of $j_n(t)$ in
   A/m^2 since $`j_{\text{ref}} = 1`$ A/m^2:
   ```math
   V(t) = V_D + \frac{j_n(t)}{j_{\text{ref}}}V_{\text{ref}},\qquad
   \nabla V(t) = \frac{j_n(t)}{j_{\text{ref}}}\nabla V_{\text{ref}}.
   ```
   Only $\nabla V$ enters the spin accumulation. The imposed potential $V_D$, a constant on a
   single surface, only shifts the potential: it is added afterwards, and the results do not
   depend on it. Solving the reference problem with $V = V_D$ instead of 0, and subtracting
   $V_D$, would not do: the reference potential (about 5e-14 V for 1 A/m^2 in a cylinder of
   200 nm with $\sigma$ = 4e6 S/m) would be drowned in the error of the solve around $V_D$ (about
   3e-8 V for $V_D$ = 1 V), giving tens of kV instead of 0.05 V once multiplied by 1e12.

   A numeric `jn` follows the same path, with the factor 1: the potential is solved with the
   `jn` of the settings and V = 0, and $V_D$ is added, $V_0 = V_D + V$. Solving around
   $V_D \neq 0$ also degraded it, less severely: with $V_D$ = 1 V on the domain wall case, the
   error of about 1e-8 V on $V - V_D$ gave 6.6e-4 (relative) on the spin accumulation.

3. **Quasi-static spin accumulation.** The spin accumulation is solved at each time step, with
   the current magnetization and the current of that time: its diffusion and relaxation times are
   assumed much shorter than the time scales of the current and of the magnetization. In the weak
   form (`spinAcc::solve`), the right-hand side of the equation of the test function $a_i$ has
   two source terms: the current polarized by the magnetization,
   ```math
   \int_\Omega \frac{\mu_B}{e}P\sigma(\nabla V\cdot\nabla a_i)\mathbf{u},
   ```
   and the spin injection through the surface $\Gamma$ with `jn` and `uP`,
   ```math
   -\int_\Gamma Q_n a_i,\qquad Q_n = -\frac{\mu_B}{e}j_n\mathbf{u}_P.
   ```
   Both are linear in $(\nabla V, j_n)$: they are computed with $`\nabla V_{\text{ref}}`$ and
   $`j_{\text{ref}} = 1`$ A/m^2, multiplied by the dimensionless factor $`j_n(t)/j_{\text{ref}}`$
   (`spinAcc::jn_factor`; 1 for a numeric `jn`, whose results are thus unchanged).

4. **Imposed spin accumulation unchanged.** The value `s` of the surface regions (Dirichlet
   condition, a reservoir) does not follow the current. If the surface with `jn` also has `s`, the
   Dirichlet condition holds there and the injection through it has no effect. With
   $\mathbf{s} = 0$ imposed and $j_n(t) = 0$ (e.g. at the start of a ramp), the right-hand side
   vanishes: the solver converges at once (zero residual) with $\mathbf{s} = 0$.

5. **Time discretization.** The spin accumulation is solved after each time step, with the new
   magnetization $`\mathbf{u}_{n+1}`$ and $`j_n(t_{n+1})`$, and is used by the next step, which
   starts at $`t_{n+1}`$ and evaluates the applied field at $`t_{n+1}`$: the current density and
   the field are evaluated at the same time, the start of the step (st-feeLLGood evaluates both at
   $t_n + \theta dt$). The time step is not adapted to the current: give
   `time_integration.max(dt)` well below its rise and fall times, or its period for an
   oscillating current, otherwise a fast variation is seen with a delay of one step, and a short
   pulse can be missed altogether.

6. **Physical time.** $t$ is the physical time in seconds, the one of `Bext` and of the `.evol`
   file. On a restart from a `.sol` file or with `initial_time`, the current goes on from that
   time.

## Implementation

- `Triangle::prm::jn_current_density` marks the surfaces with `jn: current_density`, whose `jn` is
  then the reference 1 A/m^2.
- `ExpressionParser::define_global_function` defines `current_density` as a global JavaScript
  function in the interpreters of the time functions of the field (`field_parser`,
  `field_time_parser`): each `ExpressionParser` has its own interpreter, the function has to be
  defined in each of them. As JavaScript resolves global names at the call, the order of the keys
  `current_density` and `Bext` in the settings does not matter.
- `Settings::getCurrentDensity(t)` evaluates $j_n(t)$, `Settings::hasCurrentDensity()`,
  `Settings::usesCurrentDensity()` (a surface has `jn: current_density`), and
  `Settings::jnFactor(t)`, the dimensionless factor $`j_n(t)/j_{\text{ref}}`$ if a surface has
  `jn: current_density`, 1 otherwise.
- `electrostatSolver::solve` always solves with V = 0 on the surface with `V` (reference problem
  for `jn: current_density`); `spinAcc::potential()` adds $V_D$, as `electrostatSolver::save`
  for a numeric `jn` ($V_0$ in `_V.sol`).
- `Fem::compute_all` sets `spinAcc::jn_factor` before every solve of the spin accumulation; it is
  called after the time increment, at the time of the new magnetization.

## Validation

Domain wall in a cylinder (D = 50 nm, L = 1000 nm, spin accumulation enabled, 1e-12 s), numeric
`jn` compared with `jn: current_density` and a constant `current_density`:

| Comparison | Max difference |
|---|---|
| `jn: 1e12` / `current_density: 1e12` | <Mx> 1.1e-16, <My> 8.9e-17, E_tot 6.2e-32 J |
| `jn: 5e11` / `current_density: 5e11` | <Mx> 1.1e-16, <My> 4.5e-17, E_tot 6.8e-32 J |

The results are not bitwise identical, since the potential is solved with jn = 1 A/m^2 and
multiplied, instead of being solved with the value of `jn`: the relative difference on
<My> (about 1.2e-6) is about 7e-11, far below the tolerance of the solvers. A numeric `jn` goes
through a factor 1: its results are unchanged.

With a linear rise of the current density over 5e-13 s,
`current_density: 1e12*min(t/5e-13, 1)`, and `Bext: ["1e-14*current_density(t)", 0, 0]`, the
column `current_density` follows the rise and the Zeeman energy follows the field.

Ellipsoid of `examples/ellipsoid.msh`, **spin accumulation disabled**,
`current_density: 1e12*min(t/5e-12, 1)` and the applied field `1e-13*current_density(t)` along x
(0.1 T on the plateau), given in the three forms of `Bext`: the three runs give the same <Mx>,
<My> and E_zeeman to 1e-9 (relative), `Hx` = 79577 A/m on the plateau. Without the key
`current_density`, the run stops with the `ReferenceError` above.

The two checks: `jn: current_density` without the key stops the run with the error above; the
key with a numeric `jn` gives the warning and the run goes on.

Several surfaces: cylinder of `simulations/cyl_oersted_current_density_ramp1ns` with the spin
accumulation, the current injected through the bottom face and the lateral surface: both with
`jn: current_density` and `current_density: -1e12`, or both with `jn: -1e12`, give the same
potential (to the 7 digits of the `.sol`, 0.2532 V on the bottom face, the 1D estimate
$`(j/\sigma)(L + L^2/R)`$ giving 0.25 V), the same spin accumulation (1e-6, relative) and the same
magnetization after one step (1e-12); mixing them stops the run with the error above.
