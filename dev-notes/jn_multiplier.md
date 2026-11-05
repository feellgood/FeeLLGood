# Time-dependent normal current densities: `jn_multiplier`

The key `spin_accumulation.jn_multiplier` gives a dimensionless function of time, $`f(t)`$, by
which all the fixed normal current densities `jn` of the surface regions are multiplied: the
normal current density of a surface with `jn` is $`f(t) j_n`$ at time $`t`$, in seconds. It
defaults to the constant value 1. The expressions of `Bext` can call it, e.g. for the Oersted
field of the current. In st-feeLLGood, the same current is `pulse(t)` × `J_n`, `pulse(t)` being
defined in the case library (`pulse.cc`).

## Usage

### Definition

`jn_multiplier` is either an expression of `t`, which feeLLGood wraps into
`function(t) { return (...); }`, or a JavaScript function expression of the single parameter `t`:
```yaml
mesh:
  surface_regions:
    inlet:
      jn: -1e12            # A/m^2, outward normal
      uP: [0, 1, 0]
    outlet:
      V: 0
      s: [0, 0, 0]

spin_accumulation:
  enable: true
  jn_multiplier: "min(t/1e-10, 1)"   # linear rise over 0.1 ns, then constant
```
The functions of the Math object are available without the `Math.` prefix (`min`, `exp`, `sin`,
`PI`...), as in `Bext`. The shape is arbitrary, and so is the sign. A trapezoid, the multi-line
form of YAML keeping the function readable:
```yaml
spin_accumulation:
  enable: true
  jn_multiplier: |
    function(t) {
      const t_on = 1e-10, rise = 5e-11, length = 1e-9, fall = 5e-11;  // [s]
      if (t < t_on) return 0;
      if (t < t_on + rise) return (t - t_on)/rise;
      if (t < t_on + rise + length) return 1;
      if (t < t_on + rise + length + fall) return 1 - (t - t_on - rise - length)/fall;
      return 0;
    }
```
and a sinusoidal current of 1 GHz switched on with a linear rise over 0.5 ns:
`jn_multiplier: "min(t/5e-10, 1)*sin(2*PI*1e9*t)"`.

Several surfaces may have `jn`, each with its own value: all of them are multiplied by the same
$`f(t)`$. The electrostatic problem requires at least one surface with `jn` and a single surface
with `V`.

### Sign convention

`jn` is $`\mathbf{j}\cdot\mathbf{n}`$ with $`\mathbf{n}`$ the **outward** normal of the surface: a
positive value makes the current leave the sample through that surface. To inject a current
through the face $`z = z_{\min}`$ of a wire along $`z`$ (outward normal $`-\mathbf{e}_z`$), so
that it flows along $`+z`$, `jn` is negative, and the current density in the wire is
$`j_z(t) = -f(t) j_n`$.

### Applied field following the current

The three forms of `Bext` can call `jn_multiplier(t)`, whether the spin accumulation is enabled
or not:
```yaml
# function of t
Bext: "function(t) { return [0.1*jn_multiplier(t), 0, 0]; }"

# three expressions of t
Bext: ["0.1*jn_multiplier(t)", 0, 0]

# product of a space function and a time function: only the time function can call jn_multiplier
Bext:
  space: "function(x, y, z) { return [-y, x, 0]; }"
  time: "function(t) { return 2*PI*1e-7*1e12*jn_multiplier(t); }"
```
The last form is the natural one for the Oersted field: its space profile is computed once on the
Gauss points, and its amplitude follows the current. In a wire along $`z`$ carrying the uniform
current density $`j_z`$, the Oersted field is $`\mathbf{B} = \frac{\mu_0 j_z}{2}(-y, x, 0)`$, in
T with $`x`$, $`y`$ in m. With `jn: -1e12` on the face $`z = z_{\min}`$,
$`j_z(t) = 10^{12} f(t)`$ A/m^2 and $`\mu_0/2 = 2\pi \times 10^{-7}`$, hence the time function
above: the value of `jn` has to be repeated there. feeLLGood does not compute the Oersted field of
the current itself: the profile is the user's responsibility.

If `Bext` calls `jn_multiplier` while the key is not defined, the run stops at the first
evaluation of the field with the message
`Script ReferenceError: identifier 'jn_multiplier' undefined`.

### Outputs

- `.evol`: column `jn_multiplier`, to be added to `outputs.evol_columns`, value of $`f(t)`$ at the
  time of the line (1 if the key is not defined). The columns `Hx`, `Hy`, `Hz` give the applied
  field, current included (form `Bext` function of `t` or three expressions only).
- `.sol` (with the spin accumulation): the columns are `idx mx my mz phi sx sy sz V`, the last one
  being the electrostatic potential at the time of the file, $`V(t) = V_D + f(t) V_1`$ (see the
  hypotheses). It costs no solve. When a `.sol` file is used as the initial magnetization
  (restart), the columns after `phi` are skipped: the spin accumulation and the potential are
  computed again at the start.
- `spin_accumulation.V_file`: the potential at the start, ignoring `jn_multiplier` (i.e. for
  $`f = 1`$), $`V_D + V_1`$, written once in `<file_basename>_V.sol` of `outputs.directory`.

### From Python

`settingsMaker` has no dedicated method: set the key of the dictionary directly,
```python
mySettings["spin_accumulation"]["jn_multiplier"] = "min(t/1e-10, 1)"
```

## Hypotheses

### 1. A single time profile

All the surfaces with `jn` follow the same $`f(t)`$. Surfaces with
independent profiles (e.g. a constant current on one surface and a pulse on another) would
need one solve of the potential per profile (superposition) and are not possible.

### 2. Linear and quasi-static electrokinetics

The potential is the solution of
```math
\nabla\cdot(\sigma\nabla V) = 0,\qquad
-\sigma\frac{\partial V}{\partial n} = f(t)\,j_n \text{ on the surfaces with } j_n,\qquad
V = V_D \text{ on the surface with } V.
```
The conductivity $`\sigma`$ does not depend on the magnetization (no anisotropic
magnetoresistance), on the current or on the temperature (no Joule heating), and the current
follows $`f(t)`$ instantaneously (no capacitive or inductive effect). The problem being linear,
```math
V(t) = V_D + f(t)\,V_1,\qquad \nabla V(t) = f(t)\,\nabla V_1,
```
where $`V_1`$ is solved once, at the start, with the `jn` of the settings and **V = 0** on the
surface with `V`. The imposed potential $`V_D`$, a constant on a single surface, only shifts
the potential: it is added afterwards. Solving around $`V_D \neq 0`$ would add an error of the
order of the tolerance of the solver times $`V_D`$: with $`V_D`$ = 1 V on the domain wall case
below, the error of about 1e-8 V gave 6.6e-4 (relative) on the spin accumulation.

### 3. Quasi-static spin accumulation

The spin accumulation is solved at each time step, with
the magnetization and the current of that step: its diffusion and relaxation times are assumed
much shorter than the time scales of the current and of the magnetization. In the weak form
(`spinAcc::solve`), the right-hand side of the equation of the test function $`a_i`$ has two
source terms, the current polarized by the magnetization and the spin injection through the
surfaces $`\Gamma`$ with `jn` and `uP`:
```math
\int_\Omega \frac{\mu_B}{e}P\sigma(\nabla V\cdot\nabla a_i)\,\mathbf{u}
\qquad\text{and}\qquad
-\int_\Gamma Q_n a_i,\quad Q_n = -\frac{\mu_B}{e}f(t)\,j_n\mathbf{u}_P.
```
Both are linear in $`f(t)`$: they are computed with $`\nabla V_1`$ and $`j_n`$, and multiplied
by $`f(t)`$ (`spinAcc::jn_multiplier`).

### 4. Imposed spin accumulation unchanged

The value `s` of the surface regions (Dirichlet
condition, a reservoir) does not follow the current. With $`\mathbf{s} = 0`$ imposed and
$`f(t) = 0`$ (e.g. at the start of a ramp), the right-hand side vanishes: the solver converges
at once with $`\mathbf{s} = 0`$.

### 5. Time discretization

The spin accumulation is solved after each time step, with the new
magnetization and $`f(t_{n+1})`$, and is used by the next step, which evaluates the applied
field at its start, $`t_{n+1}`$ (st-feeLLGood evaluates both at $`t_n + \theta dt`$). The time
step is not adapted to the current: give `time_integration.max(dt)` well below its rise and
fall times, or its period, otherwise a short pulse can be missed altogether.

### 6. Physical time

$`t`$ is the physical time in seconds, the one of `Bext` and of the `.evol`
file. On a restart from a `.sol` file or with `initial_time`, the current goes on from that
time.

## Implementation

- `Settings::jnMultiplier(t)` returns $`f(t)`$ if `spin_accumulation.jn_multiplier` is defined, 1
  otherwise; `ExpressionParser::is_set()` tells whether its parser has received a function.
- `ExpressionParser::define_global_function` defines `jn_multiplier` as a global JavaScript
  function in the interpreters of the time functions of the field (`field_parser`,
  `field_time_parser`): each `ExpressionParser` has its own interpreter. As JavaScript resolves
  global names at the call, the order of the keys `Bext` and `spin_accumulation` does not matter.
- `electrostatSolver::solve` solves with V = 0 on the surface with `V`; `spinAcc::potential()`
  and `electrostatSolver::save` add $`V_D`$.
- `Fem::compute_all` sets `spinAcc::jn_multiplier` before every solve of the spin accumulation;
  it is called after the time increment, at the time of the new magnetization.

## Validation

Domain wall in a cylinder (D = 50 nm, L = 1000 nm, spin accumulation enabled, 3e-12 s), compared
with the former implementation (`jn: 1e12`, key `current_density`, commit ddff946):

| Run | Max difference |
|---|---|
| `jn: 1e12`, no `jn_multiplier` | <Mx> 1.9e-16 (run-to-run noise), `.sol` identical |
| `jn: 1`, `jn_multiplier: "1e12"` | <Mx> 5.8e-16, s 2.2e-7 (relative), V identical |
| `jn: 2e12`, `jn_multiplier: 0.5` | <Mx> 1.9e-16, s 6.5e-12 (relative), V identical |
| `jn: 1e12`, `jn_multiplier: "min(t/5e-13, 1)"`, `Bext: ["1e-2*jn_multiplier(t)", 0, 0]`, compared with `jn: current_density`, `current_density: "1e12*min(t/5e-13, 1)"`, `Bext: ["1e-14*current_density(t)", 0, 0]` | <Mx> 5.0e-16, s 2.2e-7 (relative), V identical, same `Hx` |

The differences on s come from the solve of the potential with a different `jn` (the tolerance
of the solver being relative). With `spin_accumulation.V_file: true` and
`outputs.directory: out`, the file is written in `out/`.

Ellipsoid of `examples/ellipsoid.msh`, **spin accumulation disabled**,
`jn_multiplier: "min(t/5e-12, 1)"` and the applied field `0.1*jn_multiplier(t)` T along x, given
in the three forms of `Bext`: the three runs give the same <Mx>, <My> and E_zeeman to 1e-9
(relative), `Hx` = 79577 A/m on the plateau (NaN for the form space × time, as without
`jn_multiplier`).
