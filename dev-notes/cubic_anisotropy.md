# Cubic anisotropy in feeLLGood: derivative of the effective field with respect to the magnetization

## Context

LLG is solved with the θ-scheme of
F. Alouges, E. Kritsikis, J. Steiner, J.-C. Toussaint,
_A convergent and precise finite element scheme for Landau–Lifschitz–Gilbert equation_,
Numer. Math. **128**(3), 407–430 (2014).
At each step, one looks for $\mathbf{v}$ in the plane tangent to $\mathbf{u}^n$ such that, for every
tangent test function $\mathbf{w}$,
```math
\int (\tilde\alpha\,\mathbf{v} + \mathbf{u}^n\times\mathbf{v})\cdot\mathbf{w}
 + \theta\,\delta t\,(1+R)\,a(\mathbf{v},\mathbf{w})
 = \int \mathbf{H}_{\rm eff}(\mathbf{u}^n)\cdot\mathbf{w}
 + \theta\,\delta t\int \mathbf{H}'_{\rm eff}(\mathbf{u}^n)[\mathbf{v}^{n-1}]\cdot\mathbf{w},
\qquad (1)
```
(up to the factor $\gamma_0$ between time and reduced time),
then $\mathbf{u}^{n+1}=(\mathbf{u}^n+\delta t\mathbf{v})/|\mathbf{u}^n+\delta t\mathbf{v}|$.
The last term is a Taylor expansion of the field at $t_n+\theta\delta t$:
```math
\mathbf{H}_{\rm eff}(\mathbf{u}^n+\theta\,\delta t\,\mathbf{v}) =
   \mathbf{H}_{\rm eff}(\mathbf{u}^n) + \theta\,\delta t\,\mathbf{H}'_{\rm
 eff}(\mathbf{u}^n)[\mathbf{v}] + O(\delta t^2),
```
where
```math
\mathbf{H}'_{\rm eff}(\mathbf{u})[\mathbf{v}] =
   \frac{\partial \mathbf{H}_{\rm eff}}{\partial \mathbf{u}}\,\mathbf{v}
```
is the derivative of the field with respect to the magnetization in the direction
$\mathbf{v}$. It is this term that makes the scheme of second order in time; it must be the _exact_
derivative of the field. In the code, the anisotropy fields are computed in `Tet::calc_aniso_uniax`
and `Tet::calc_aniso_cub` (`src/tetra.cpp`), called from `Tet::integrales`: they append to
`H_aniso`, at each Gauss point, the field $\mathbf{H}(\mathbf{u})$ plus the second order term
`s_dt`$\mathbf{H}'(\mathbf{u})[\mathbf{v}]$, with `s_dt` $=\theta\delta t$ (`THETA` is defined
in `config.h`).

## Cubic anisotropy: energy and field

Let $\mathbf{k}_1,\mathbf{k}_2,\mathbf{k}_3$ be the orthonormal cubic axes (`ex`, `ey`, `ez` of the
material in the settings file, `Tetra::prm::ex,ey,ez` in the code) and
$a_l=\mathbf{k}_l\cdot\mathbf{u}$ the direction cosines, with $\sum_l a_l^2=|\mathbf{u}|^2=1$.
The energy density is (`Tet::cubicAnisotropyEnergy`)
```math
e(\mathbf{u}) = K_3\,\bigl(a_1^2a_2^2+a_2^2a_3^2+a_3^2a_1^2\bigr).
```
Since $`\partial a_l/\partial\mathbf{u}=\mathbf{k}_l`$ and
$`\partial e/\partial a_l = 2K_3a_l\sum_{m\neq l}a_m^2 = 2K_3a_l(1-a_l^2)`$ on the unit sphere,
the field $\mathbf{H}=-\frac{1}{\mu_0 M_s}\frac{\partial e}{\partial\mathbf{u}}$ is
```math
\mathbf{H}(\mathbf{u}) = -\frac{2K_3}{\mu_0 M_s}\sum_{l=1}^{3} a_l\,(1-a_l^2)\,\mathbf{k}_l .
\qquad (2)
```
In the code, `K3bis` $=2K_3/(\mu_0 M_s)$. This expression of the field is unchanged by the
correction, as well as the contribution $\mathbf{u}\cdot\mathbf{H}$ returned by `calc_aniso_cub`
(used to compute the stabilizing effective damping).

## Exact derivative (new expression)

Differentiating (2) with
$\frac{\partial a_l}{\partial\mathbf{u}}\cdot\mathbf{v}=\mathbf{k}_l\cdot\mathbf{v}$:
```math
\boxed{\;
\mathbf{H}'(\mathbf{u})[\mathbf{v}] = -\frac{2K_3}{\mu_0
M_s}\sum_{l=1}^{3}(\mathbf{k}_l\cdot\mathbf{v})\,\bigl(1-3a_l^2\bigr)\,\mathbf{k}_l
\;}
\qquad (3)
```
that is, in matrix form,
```math
\frac{\partial\mathbf{H}}{\partial\mathbf{u}}
 = -\frac{2K_3}{\mu_0 M_s}\sum_{l=1}^{3}\bigl(1-3a_l^2\bigr)\,\mathbf{k}_l\otimes\mathbf{k}_l ,
\qquad
\frac{\partial H_d}{\partial u_{d'}}
 = -\frac{2K_3}{\mu_0 M_s}\sum_{l}\bigl(1-3a_l^2\bigr)\,k_{l,d}\,k_{l,d'} .
```
This matrix is symmetric, as it must be: up to the factor $-1/(\mu_0 M_s)$, it is the Hessian of the
energy. In the eigenbasis $(\mathbf{k}_1,\mathbf{k}_2,\mathbf{k}_3)$ it is diagonal, with
eigenvalues $-\frac{2K_3}{\mu_0 M_s}(1-3a_l^2)$.

### Remark (constraint $|\mathbf{u}|=1$)

Differentiating instead the unconstrained form $a_l\sum_{m\ne l}a_m^2$ adds the term
$-\frac{4K_3}{\mu_0 M_s}\sum_l a_l(\mathbf{u}\cdot\mathbf{v})\mathbf{k}_l$, which vanishes for
$\mathbf{v}$ tangent ($\mathbf{u}\cdot\mathbf{v}=0$), the only case used by the scheme. Both forms
therefore coincide.

## Former expression (erroneous)

The former code of `Tet::calc_aniso_cub`, inherited from the reference code
(`MuMag_integrales.cc`), was
```c++
Eigen::Vector3d tmp = uk_v.cwiseProduct(ex);
H_aniso.col(npi) += -K3bis * (uk_uuu(0) * ex + uk_uuu(1) * ey + uk_uuu(2) * ez
          + s_dt * tmp.cwiseProduct( Eigen::Vector3d(1, 1, 1)
              - 3*uk_u.cwiseProduct(uk_u) ));
```
with `uk_u(l)` $=a_{l+1}$ and `uk_v(l)` $`=\mathbf{k}_{l+1}\cdot\mathbf{v}`$. In the reference
code, the same term read `Ht[d] = -2*K3/Js* uk_dv*(1-3*uk_du*uk_du)*uk0d`, where
`uk0d` $`=k_{1,d+1}`$ is component $d$ of the _first_ axis. Written with indices $d=1,2,3$:
```math
H'_{\rm old,\,d}(\mathbf{u})[\mathbf{v}] = -\frac{2K_3}{\mu_0
M_s}\,(\mathbf{k}_d\cdot\mathbf{v})\,\bigl(1-3a_d^2\bigr)\,k_{1,d},
\qquad
\frac{\partial\mathbf{H}_{\rm old}}{\partial\mathbf{u}}
 = -\frac{2K_3}{\mu_0 M_s}
   \sum_{d=1}^{3}\bigl(1-3a_d^2\bigr)\,k_{1,d}\;\mathbf{e}_d\otimes\mathbf{k}_d .
\qquad (4)
```
Compared with (3), two errors combine:
1. the index of the axis $l$ is identified with the Cartesian component $d$: each component only
   keeps one axis, instead of the sum over the three axes;
2. the component taken is the one of the first axis, $k_{1,d}$, instead of $k_{d,d}$ (the factor
   `ex`, or `uk00`, `uk01`, `uk02` in the reference code, was copied from the uniaxial term, whose
   axis is $\mathbf{k}_1$).
Even for axes aligned with the Cartesian frame ($`\mathbf{k}_l=\mathbf{e}_l`$), where the correct
derivative is $`H'_d=-\frac{2K_3}{\mu_0 M_s}v_d(1-3u_d^2)`$, the former expression only keeps the
component $d=1$ ($`k_{1,d}=\delta_{1d}`$): the second-order correction of components 2 and 3 is
missing. For general axes, the matrix (4) is not even symmetric.

## Consequences and correction

The error only affects the term $\theta\delta t\mathbf{H}'[\mathbf{v}]$ of (1):
the scheme remains consistent and stable, but the local error is $O(\delta t^2)$
instead of $O(\delta t^3)$ for the cubic contribution, so the scheme loses its second
order in time when $K_3\neq0$. Nothing changes when $K_3=0$.

The correction (`Tet::calc_aniso_cub`, `src/tetra.cpp`) implements (3), with the cubic axes as the
columns of a matrix $k=[\mathbf{k}_1\ \mathbf{k}_2\ \mathbf{k}_3]$:
```c++
Eigen::Matrix3d k;
k << ex, ey, ez;
...
Eigen::Vector3d uk_u = k.transpose() * U.col(npi);
Eigen::Vector3d uk_v = k.transpose() * V.col(npi);
Eigen::Vector3d uk_uuu = uk_u.unaryExpr( [](const double x){ return x*(1.0 - x*x);} );
Eigen::Vector3d dH = uk_v.cwiseProduct( Eigen::Vector3d(1, 1, 1)
                                        - 3*uk_u.cwiseProduct(uk_u) );
H_aniso.col(npi) += -K3bis * k * (uk_uuu + s_dt * dH);
```

The unit tests (`unit-tests/ut_anisotropy.cpp`) check it in two ways:
* `anisotropy_cubic` compares `H_aniso` with an explicit expression
  of $\mathbf{H}+\theta\delta t\mathbf{H}'[\mathbf{v}]$, where the former reference expression
  of $\mathbf{H}'$ has been replaced by (3);
* `anisotropy_cubic_derivative` checks
  $\mathbf{H}'[\mathbf{v}]=\partial\mathbf{H}/\partial\mathbf{u}\cdot\mathbf{v}$ against a centered
  finite difference of the field, for random axes (relative deviation $<10^{-6}$; about $10^{-10}$
  in practice). With the former expression, the relative deviation was of order $1$.

## Effect on a simulation

The following figures were obtained with st-feeLLGood, which had the same erroneous
expression and received the same correction; they have not been reproduced with
feeLLGood. Regression case `short_cyl_rich` (cylinder, $K_3=5\times10^4 \mathrm{J}/\mathrm{m}^3$,
together with uniaxial anisotropy, DMI, applied field and spin transport; final time 4.52 ps).
Maximum deviation of the average magnetization $\langle\mathbf{u}\rangle$ over the 11 outputs:

| comparison                       | $\max\lvert\Delta\langle\mathbf{u}\rangle\rvert$ |
|----------------------------------|--------------------------------------------------|
| former vs corrected $\mathbf{H}'$ (default time step)          | $4.6\times10^{-7}$ |
| former vs corrected $\mathbf{H}'$ (step ÷4, max(_δu_) = 0.005) | $1.6\times10^{-7}$ |
| default step vs step ÷4 (former and corrected)                 | $1.3\times10^{-4}$ |

On this case, the correction changes the result by much less than the total time
discretization error, which is dominated by the other terms (the cubic anisotropy is weak
compared with the demagnetizing field): the effect of the correction is not visible on
the global convergence. It becomes significant when the cubic anisotropy dominates the
effective field (large $K_3$, weak shape anisotropy).
