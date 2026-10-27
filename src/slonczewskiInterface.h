#ifndef SLONCZEWSKI_INTERFACE_H
#define SLONCZEWSKI_INTERFACE_H

/** \file slonczewskiInterface.h
\brief effective boundary condition of the spin diffusion problem for a spin injection through a
polarizer F1 / thin non-magnetic spacer N / free layer F2 stack, F1 and N being not meshed (option 2
of the Slonczewski note, magnetoelectronic circuit theory).

The interface N/F with the magnetization \f$ \mathbf{u} \f$ is described by the 3x3 conductance
matrix \f$ G(\mathbf{u}) = G_\parallel\,\mathbf{u}\mathbf{u}^T + G_r(I - \mathbf{u}\mathbf{u}^T)
- G_i\,[\mathbf{u}]_\times \f$, and its polarization \f$ P \f$. The spin current entering F2 is
\f$ \mathbf{q} = K\,P_{\mathrm{eff}}(\mathbf{m}_0) + G_s(\mathbf{s}_{F1} - \mathbf{s}_F) \f$, with
\f$ K = (\mu_B/e) J \f$, \f$ J > 0 \f$ for electrons flowing from F1 to F2, \f$ \mathbf{m}_0 \f$
the magnetization of F2 on the interface and \f$ \mathbf{s}_F \f$ its spin accumulation. The bulk
of F1 is folded into the F1/N interface (renormalizeF1()), so that \f$ \mathbf{s}_{F1} = 0 \f$.

The settings (StackSettings) give the stack in physical units: resistance-area products AR
(Ohm m^2) and asymmetries gamma of the interfaces (Valet-Fert), mixing conductances G_up,down of
Brataas et al. (S/m^2), conductivity of F1. makeStack() converts them into the conductances of the
spin diffusion problem (m/s, Stack) with the factor k = 2/(e^2 N0) of spinFactor(), N0 being the
density of states of the sample on the surface. With D0 = 2 sigma/(e^2 N0), the spin accumulation of
feeLLGood is s = mu_B N0 mu_s/2, mu_s = (mu+ - mu-)/2 being the spin accumulation of Valet-Fert and
of Brataas et al. (validation/valet_fert_1d). The spin current through an interface being
(1 - gamma^2)(G_up + G_down) mu_s = mu_s/AR along u and 2 G_up,down mu_s across u (circuit theory),
Gpar = k/AR and Gr + i Gi = 2 k G_up,down: for a non-magnetic interface (gamma = 0,
G_up = G_down = G_up,down), the conductance is isotropic.
*/

#include <cmath>

#include <eigen3/Eigen/Dense>

#include "config.h"

namespace Slonczewski
    {
/** conductance matrix of an interface N/F, the ferromagnet having the magnetization u (unit
 * vector): G u = Gpar u, and for t orthogonal to u, G t = Gr t - Gi u x t */
inline Eigen::Matrix3d conductance(const Eigen::Vector3d &u /**< [in] magnetization */,
                                   const double Gpar /**< [in] longitudinal conductance */,
                                   const double Gr /**< [in] real part of the mixing conductance */,
                                   const double Gi /**< [in] imaginary part of the mixing conductance */)
    {
    Eigen::Matrix3d u_cross;  // u_cross w = u x w
    u_cross << 0, -u.z(), u.y(),
               u.z(), 0, -u.x(),
               -u.y(), u.x(), 0;
    const Eigen::Matrix3d uuT = u * u.transpose();
    return Gpar * uuT + Gr * (Eigen::Matrix3d::Identity() - uuT) - Gi * u_cross;
    }

/** factor converting a conductance per unit area (S/m^2), relating the spin current to mu_s, into a
 * conductance of the spin diffusion problem (m/s), s = mu_B N0 mu_s/2: 2/(e^2 N0) */
inline double spinFactor(const double N0 /**< [in] density of states of the sample, J^-1 m^-3 */)
    { return 2.0 / (CHARGE_ELECTRON * CHARGE_ELECTRON * N0); }

/** longitudinal conductance per unit area (S/m^2) of the bulk of F1, of thickness t, with a
 * reservoir (s = 0) on its far side: sigma (1 - P^2)/l_sf coth(t/l_sf) (beta = beta' = P, l_sf
 * being the longitudinal spin diffusion length); t = infinity is allowed */
inline double bulkConductance(const double sigma /**< [in] conductivity of F1, S/m */,
                              const double P /**< [in] bulk polarization rate of F1 */,
                              const double lsf /**< [in] spin diffusion length of F1 */,
                              const double t /**< [in] thickness of F1 */)
    { return sigma * (1.0 - P * P) / (lsf * std::tanh(t / lsf)); }

/** \struct F1prm
polarization and longitudinal conductance of the F1/N interface */
struct F1prm
    {
    double P;    /**< polarization */
    double Gpar; /**< longitudinal conductance */
    };

/** renormalization of the F1/N interface by the bulk of F1, its spin accumulation being collinear
 * to p: the spin current of the interface keeps the form of s_F1 = 0 with
 * P* = (beta Gpar + P gF)/(gF + Gpar) and Gpar* = gF Gpar/(gF + Gpar) (in series); the mixing
 * conductances are unchanged */
inline F1prm renormalizeF1(const double beta /**< [in] bulk polarization rate of F1 */,
                           const double gF /**< [in] bulk conductance of F1, bulkConductance() */,
                           const F1prm &itf /**< [in] F1/N interface */)
    {
    const double sum = gF + itf.Gpar;
    return {(beta * itf.Gpar + itf.P * gF) / sum, gF * itf.Gpar / sum};
    }

/** \struct BC
effective boundary condition: entering spin current q = K Peff + Gs (s_F1 - s_F) */
struct BC
    {
    Eigen::Vector3d Peff; /**< effective polarization vector, dimensionless */
    Eigen::Matrix3d Gs;   /**< conductance of the two interfaces in series */
    };

/** effective boundary condition, the spin accumulation s_N of the spacer being eliminated from the
 * conservation of the spin current through N: with A = G1 + G2,
 * Peff = G2 A^-1 P1 p + G1 A^-1 P2 m0 and Gs = G2 A^-1 G1 */
inline BC effectiveBC(const Eigen::Matrix3d &G1 /**< [in] F1/N conductance, G(p) */,
                      const double P1 /**< [in] F1/N polarization */,
                      const Eigen::Vector3d &p /**< [in] magnetization of F1 */,
                      const Eigen::Matrix3d &G2 /**< [in] N/F2 conductance, G(m0) */,
                      const double P2 /**< [in] N/F2 polarization */,
                      const Eigen::Vector3d &m0 /**< [in] magnetization of F2 on the interface */)
    {
    const Eigen::PartialPivLU<Eigen::Matrix3d> A(G1 + G2);
    return {G2 * A.solve(P1 * p) + G1 * A.solve(P2 * m0), G2 * A.solve(G1)};
    }

/** \struct Stack
parameters of the polarizer stack of a surface: the F1/N interface (renormalized by the bulk of F1,
see renormalizeF1()) and, optionally, the N/F2 interface. Without N/F2 (spacer meshed, or perfect
N/F2 interface), the condition does not depend on the magnetization of the sample: Peff = P1 p and
Gs = G1 */
struct Stack
    {
    Eigen::Vector3d p = Eigen::Vector3d(0, 0, 1); /**< magnetization of F1, unit vector */
    double P1 = 0;     /**< polarization of F1/N */
    double G1par = 0;  /**< longitudinal conductance of F1/N */
    double G1r = 0;    /**< real part of the mixing conductance of F1/N */
    double G1i = 0;    /**< imaginary part of the mixing conductance of F1/N */
    bool hasNF2 = false; /**< true if the N/F2 interface is given */
    double P2 = 0;     /**< polarization of N/F2 */
    double G2par = 0;  /**< longitudinal conductance of N/F2 */
    double G2r = 0;    /**< real part of the mixing conductance of N/F2 */
    double G2i = 0;    /**< imaginary part of the mixing conductance of N/F2 */

    /** effective boundary condition for the magnetization m0 of the sample on the interface (unit
     * vector, not used without N/F2) */
    BC bc(const Eigen::Vector3d &m0) const
        {
        const Eigen::Matrix3d G1 = conductance(p, G1par, G1r, G1i);
        if (!hasNF2)
            { return {P1 * p, G1}; }
        return effectiveBC(G1, P1, p, conductance(m0, G2par, G2r, G2i), P2, m0);
        }
    };

/** \struct Interface
interface N/F in physical units */
struct Interface
    {
    double AR = 0;    /**< resistance-area product, Ohm m^2 */
    double gamma = 0; /**< asymmetry of the interface resistance (polarization), |gamma| <= 1 */
    double Gr = 0;    /**< real part of the mixing conductance G_up,down of Brataas et al., S/m^2 */
    double Gi = 0;    /**< imaginary part of the mixing conductance G_up,down, S/m^2 (its sign is not
                         checked by the validations) */
    };

/** \struct StackSettings
polarizer stack of a surface, in physical units, as given by the settings */
struct StackSettings
    {
    Eigen::Vector3d p = Eigen::Vector3d(0, 0, 1); /**< magnetization of F1, unit vector */
    bool hasF1 = false;       /**< true if the bulk of F1 is given */
    double sigma1 = 0;        /**< conductivity of F1, S/m */
    double P1bulk = 0;        /**< bulk polarization rate of F1 */
    double lsf1 = 0;          /**< spin diffusion length of F1, m */
    double t1 = INFINITY;     /**< thickness of F1, m */
    Interface F1N;            /**< interface F1/N */
    bool hasNF2 = false;      /**< true if the interface N/F2 is given */
    Interface NF2;            /**< interface N/F2 */
    };

/** stack in the units of the spin diffusion problem (m/s), N0 being the density of states of the
 * sample on the surface, k = 2/(e^2 N0): the longitudinal conductance of an interface is k/AR, its
 * mixing conductance 2 k G_up,down; the bulk of F1 is folded into the F1/N interface */
inline Stack makeStack(const StackSettings &st /**< [in] */,
                       const double N0 /**< [in] density of states of the sample, J^-1 m^-3 */)
    {
    const double k = spinFactor(N0);
    Stack s;
    s.p = st.p;
    s.P1 = st.F1N.gamma;
    s.G1par = k / st.F1N.AR;
    s.G1r = 2 * k * st.F1N.Gr;
    s.G1i = 2 * k * st.F1N.Gi;
    if (st.hasF1)
        {
        const double gF = k * bulkConductance(st.sigma1, st.P1bulk, st.lsf1, st.t1);
        const F1prm r = renormalizeF1(st.P1bulk, gF, {s.P1, s.G1par});
        s.P1 = r.P;
        s.G1par = r.Gpar;
        }
    s.hasNF2 = st.hasNF2;
    if (st.hasNF2)
        {
        s.P2 = st.NF2.gamma;
        s.G2par = k / st.NF2.AR;
        s.G2r = 2 * k * st.NF2.Gr;
        s.G2i = 2 * k * st.NF2.Gi;
        }
    return s;
    }
    }  // namespace Slonczewski

#endif
