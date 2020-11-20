#ifndef time_integration_h
#define time_integration_h

#include <iostream>

/** all timing parameters for integrating in time LLG with adaptative time-step */

class timing
{
    public:
        inline timing():t(0),tf(0),DTMIN(1e-14),DTMAX(1e-7),TAUR(100.*DTMAX) {set_dt_init();}
    
    double t;/**< physical current time of the simulation */
    double dt;/**< time-step */
    double tf;/**< final time of the simulation */

    /** minimum step time for time integrator */
    double DTMIN; //1e-14;

/** maximum step time for time integrator */
    double DTMAX;//  1e-5 en stat ;  1e-7 en dynamique;

/** reduced \f$ \tau_r = 100 DT_{\mathrm{max}} \f$ */
    double TAUR;
    
    /** some informations */
    inline void infos()
    {
    std::cout << "\t final time of the simulation\t\t" << tf << std::endl;
    std::cout << "\t initial time step\t" << dt << std::endl << std::endl;
    }

    /** initialize variable time step with geometric average of DTMIN and DTMAX */
    inline void set_dt_init() {dt = sqrt(DTMIN*DTMAX);}
    };



#endif
