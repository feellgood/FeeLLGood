# Spin injection through a polarizer stack: 1D validation against Valet–Fert

`valet_fert.py` checks the Robin boundary condition of the surface key `slonczewski`
(`src/slonczewskiInterface.h`, `spinAcc::stackBoundaryCondition`) on a collinear 1D case, against
the analytic solution of Valet and Fert (Phys. Rev. B **48**, 7099, 1993).

The stack is F1 (bulk, thickness $t_1$, reservoir on its far side) / interface F1/N / spacer N of
zero thickness / interface N/F2 / F2, F2 being the meshed wire (`wire.geo`, 10 × 10 × 300 nm,
h = 1 nm). The magnetizations are along $z$, F1 parallel (P) or antiparallel (AP) to F2.

It checks the conversion, done by feeLLGood (`Slonczewski::makeStack`), of the interface
resistance-area product $AR^*$ (Ω·m²) of the settings into the longitudinal conductance of the spin
diffusion problem:
```math
G_\parallel = \frac{2}{e^2N_0\,AR^*}\quad[\mathrm{m/s}],
```
with $s = \mu_B N_0\,\mu_s/2$ (A/m), $\mu_s = (\mu_+ - \mu_-)/2$, $N_0$ being the density of states of F2 (the
volume region bordering the surface): the stack is not meshed, all its conductances are converted
with it. The asymmetry `gamma` of an interface is the one of Valet–Fert.

## Usage

```
./valet_fert.py [path to feellgood]    # default ../../src/feellgood
```
The mesh is built with `gmsh` if missing. The exit status is 0 if the maximum error on $s_z$,
relative to $\max|s_z|$, is below 1 % in both configurations.

## Results (branch slonczewski)

| mesh   | P: $s_z$(inlet) feeLLGood / Valet–Fert | AP: $s_z$(inlet) feeLLGood / Valet–Fert | max error |
|--------|-----------------------------------------|------------------------------------------|-----------|
| h = 2 nm | 45.690 / 45.858 A/m                   | −163.994 / −164.596 A/m                  | 8.2e-3    |
| h = 1 nm | 45.811 / 45.858 A/m                   | −164.427 / −164.596 A/m                  | 2.6e-3    |

The error decreases as $h^2$ (P1 elements). With the factor 1 instead of 2 in $G_\parallel$, the
error is 8.5 % (P configuration): the test fails.

## Documents

`valet_fert_1d.tex` (English) and `valet_fert_1d_french.tex` (French) explain the test in detail,
for readers new to spin transport: the assumptions and equations of the Valet–Fert model, the test
case, its analytic solution, what feeLLGood solves, and the expected results. Their circuit reading
of the solution gives the polarization of the injected current as the average of the
polarizations of the elements of the stack, weighted by their spin resistances. Build them with
```
pdflatex valet_fert_1d.tex          # twice, for the table of contents
pdflatex valet_fert_1d_french.tex
```
Their figure `valet_fert_1d_profile.pdf` is built by `plot_profile.py`, after `valet_fert.py`.
