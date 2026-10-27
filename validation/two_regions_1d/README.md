# Continuity of s between regions of different densities of states: 1D validation

feeLLGood takes the spin accumulation $s$ continuous between the meshed volume regions, as is the
spin current. With $s = \mu_B N_0\,\mu_s/2$ ($\mu_s = (\mu_+-\mu_-)/2$, `validation/valet_fert_1d`),
this is the continuity of the spin chemical potential $\mu_s$ of Valet–Fert only if the regions have
the same `dens_state`; otherwise, it imposes
```math
N_{0,A}\,\mu_{s,A} = N_{0,B}\,\mu_{s,B}
```
on the interface between the regions A and B. In verbose mode, feeLLGood warns about the
neighbouring regions of different `dens_state` (`spinAcc::checkDensityOfStates`).

`two_regions.py` checks this model on a collinear 1D case, against a reference written in the
variables of Valet–Fert with the interface rule above.

The wire (`two_regions.geo`, 4 × 4 × 300 nm, h = 1 nm) is a non-magnetic spacer N (0 < z < 40 nm,
σ = 6e7 S/m, `dens_state` 0.6e47 J⁻¹m⁻³, $l_{sf}$ = 50 nm), then the free layer F2 (σ = 4e6 S/m,
`dens_state` 1.344e47, P = 0.5, $l_{sf}$ = 10 nm), then a reservoir ($s = 0$). The spin is injected
on z = 0 through a polarizer stack without N/F2 interface (the spacer being meshed): F1 (ideal
reservoir) / F1/N ($AR^*$ = 0.5 fΩ·m², $\gamma$ = 0.75). The magnetizations are along $z$, F1
parallel (P) or antiparallel (AP) to F2.

## Usage

```
./two_regions.py [path to feellgood]    # default ../../src/feellgood
```
The mesh is built with `gmsh` if missing. The exit status is 0 if the maximum error on $s_z$,
relative to $\max|s_z|$, is below 1 % in both configurations (about 8 s).

## Results (branch slonczewski)

| | $s_z$ on N/F2, feeLLGood / model ($s$ continuous) | max error |
|----|-----------------------|---------|
| P  | −7.130 / −7.131 A/m   | 1.4e-3  |
| AP | −23.365 / −23.369 A/m | 1.4e-3  |

With the same `dens_state` in both regions, the two models (continuity of $s$ or of $\mu_s$)
coincide (1e-15). Here ($N_{0,F}/N_{0,N}$ = 2.24), they differ by 100 % of $\max|s|$, on the F2 side
of the interface: with $\mu_s$ continuous, $s$ jumps by the factor 2.24 there, and the accumulation
in F2 at the interface is about twice the one of the model of feeLLGood:

| | $s$ continuous: $s_N$ / $s_F$ at the interface | $\mu_s$ continuous: $s_N$ / $s_F$ at the interface |
|----|-------------------|-------------------|
| P  | −7.13 / −7.13 A/m   | −6.40 / −14.34 A/m  |
| AP | −23.37 / −23.37 A/m | −20.98 / −47.00 A/m |

The torque in F2 depends directly on it: to have the continuity of $\mu_s$ of Valet–Fert, give the
same `dens_state` to all the regions.
