#ifndef TUNNEL_JUNCTION_H
#define TUNNEL_JUNCTION_H

/** \file tunnelJunction.h
\brief boundary conditions of the electrostatic and spin diffusion problems for a current injected
through a tunnel barrier F1 / I / F2, the electrode F1 and the barrier I being not meshed (part
"Jonction tunnel magnétique" of dev-notes/slonczewski_polarisation.tex).

The conductance per unit area of the barrier follows the cosine law, which only assumes a
spin-conserving tunneling:
\f$ G_T(\mathbf{m}_0) = (1 + \mathcal{P}\,\mathbf{p}\cdot\mathbf{m}_0)/((1 + \mathcal{P})\,RA_P) \f$,
\f$ \mathcal{P} = \mathrm{TMR}/(2 + \mathrm{TMR}) \f$, \f$ RA_P \f$ being the resistance-area
product in the parallel configuration, \f$ \mathbf{p} \f$ the magnetization of F1 and
\f$ \mathbf{m}_0 \f$ the one of the sample F2 on the surface. The normal current density leaving the
sample through the barrier is \f$ j_n = G_T(\mathbf{m}_0)\,(V - V_1) \f$, \f$ V_1 \f$ being the
potential of F1 (\f$ j_n > 0 \f$: electrons entering the sample, as in feeLLGood): a Robin condition
for the potential.

The spin current entering the sample is \f$ \mathbf{q} = (\mu_B/e)\,j_n\,\mathbf{P}_T(\mathbf{m}_0)
\f$, with the polarization of the Jullière model extended to non-collinear magnetizations,
\f$ \mathbf{P}_T = (P_1\mathbf{p} + P_2\mathbf{m}_0)/(1 + P_1P_2\,\mathbf{p}\cdot\mathbf{m}_0) \f$,
\f$ P_1 P_2 = \mathcal{P} \f$: P1 and P2 are effective tunnel polarizations, not bulk properties of
the electrodes. The feedback of the spin accumulation of the sample on the barrier is neglected,
RA being 10^3 to 10^4 times the spin resistance of the sample: the condition is a flux.
*/

#include <cmath>

#include <eigen3/Eigen/Dense>

namespace Tunnel
    {
/** \struct Settings
tunnel barrier of a surface, in physical units, as given by the settings */
struct Settings
    {
    Eigen::Vector3d p = Eigen::Vector3d(0, 0, 1); /**< magnetization of F1, unit vector */
    double RA = 0;   /**< resistance-area product in the parallel configuration, Ohm m^2 */
    double TMR = 0;  /**< tunnel magnetoresistance (G_P - G_AP)/G_AP, 1.5 for 150 % */
    double P1 = NAN; /**< effective tunnel polarization of F1, sqrt(TMR/(2 + TMR)) by default
                        (symmetric junction), see polarizations() */
    double V1 = 0;   /**< potential of F1, V, if V1_voltage is false */
    bool V1_voltage = false; /**< true if the potential of F1 is the function voltage(t) of the
                                settings (`V: voltage`) */
    };

/** potential of F1 (V): V1, or the value of voltage(t) at the current time if V1_voltage */
inline double potentialF1(const Settings &st /**< [in] */,
                          const double voltage /**< [in] voltage(t) at the current time, V */)
    { return st.V1_voltage ? voltage : st.V1; }

/** coefficient of the cosine law, P = TMR/(2 + TMR) = (G_P - G_AP)/(G_P + G_AP) */
inline double cosineCoefficient(const double TMR /**< [in] */)
    { return TMR / (2.0 + TMR); }

/** \struct Polarizations
effective tunnel polarizations of F1 and F2, P1 P2 = cosineCoefficient(TMR) */
struct Polarizations
    {
    double P1; /**< polarization of F1 */
    double P2; /**< polarization of F2 */
    };

/** effective tunnel polarizations: P1 as given (sqrt(TMR/(2 + TMR)) if NAN, symmetric junction),
 * P2 = TMR/(2 + TMR)/P1 */
inline Polarizations polarizations(const Settings &st /**< [in] */)
    {
    const double c = cosineCoefficient(st.TMR);
    const double P1 = std::isnan(st.P1) ? std::sqrt(c) : st.P1;
    return {P1, c / P1};
    }

/** conductance per unit area of the barrier (S/m^2), cosine law between G_P = 1/RA and
 * G_AP = G_P/(1 + TMR) */
inline double conductance(const Settings &st /**< [in] */,
                          const Eigen::Vector3d &m0 /**< [in] magnetization of the sample, unit
                                                       vector */)
    {
    const double c = cosineCoefficient(st.TMR);
    return (1.0 + c * st.p.dot(m0)) / ((1.0 + c) * st.RA);
    }

/** polarization of the current injected into the sample, (P1 p + P2 m0)/(1 + P1 P2 p.m0) */
inline Eigen::Vector3d polarization(const Settings &st /**< [in] */,
                                    const Eigen::Vector3d &m0 /**< [in] magnetization of the
                                                                 sample, unit vector */)
    {
    const Polarizations P = polarizations(st);
    return (P.P1 * st.p + P.P2 * m0) / (1.0 + P.P1 * P.P2 * st.p.dot(m0));
    }
    }  // namespace Tunnel

#endif
