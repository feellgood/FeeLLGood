/** \file feellgoodSettings.h
\brief many settings to give some parameters to the solver, boundary conditions for the problem, the output file format wanted by the user
*/

#ifndef fellgoodSettings_h
#define fellgoodSettings_h


/* for statics 
const double EPSILON = 1e-40;
const double DUMAX   = 0.1;
const double DTMIN   = 1e-14;
const double DTMAX   = 1.e-5;
*/



/* for dynamics */

/** \f$ \epsilon \f$ is smallest possible value (for what value ? to check) */
const double EPSILON = 1e-40;

/** min for du step */
const double DUMIN   = 1e-9;	//1e-6

/** max for du step */
const double DUMAX   = 0.02;	// 0.02

/** minimum step time for time integrator */
const double DTMIN   = 1e-14;


/** maximum step time for time integrator */
const double DTMAX   = 1e-7;	// 1e-7

/** no idea */
const double TAUR    = 100.*DTMAX;

const double mu0 = 4.*M_PI*1e-7;/**< \f$ \mu_0 = 4 \pi 10^{-7} \f$ */
const double nu0 = 1./mu0;/**< \f$ \nu_0 = 1/\mu_0 \f$ */


#endif /* feellgoodSettings_h */
