#ifndef time_integration_h
#define time_integration_h

#include <iostream>

/** all timing parameters for integrating in time LLG with adaptative time-step */

class timing
{
    public:
        /** constructor : fully initializes the timing parameters */
        inline timing(const double _t,const double _tf,const double _dtmin,const double _dtmax):tf(_tf),DTMIN(_dtmin),DTMAX(_dtmax),TAUR(100.*DTMAX),t(_t) {set_dt_init();}
    
    /** final time of the simulation */
    const double tf;

    /** minimum step time for time integrator */
    const double DTMIN; //1e-14;

/** maximum step time for time integrator */
    const double DTMAX;//  1e-5 en stat ;  1e-7 en dynamique;

/** reduced \f$ \tau_r = 100 DT_{\mathrm{max}} \f$ */
    const double TAUR;
    
    /** this prefactor must be synchronized with dt value */
    double prefactor;
    
    /** some informations */
    inline void infos()
    {
    std::cout << "\t final time of the simulation\t\t" << tf << std::endl;
    std::cout << "\t initial time step\t" << dt << std::endl << std::endl;
    }

    /** getter for time step dt */
    inline double get_dt() const {return dt;}

    /** setter for time step dt : prefactor is computed from the new dt value to maintain synchronization of dt and prefactor */
    inline void set_dt(double _dt) {dt=_dt;set_prefactor();}

    /** basic test to know if dt is too small */
    inline bool is_dt_TooSmall() const { return(dt < DTMIN);}

    /** increment time t with time step dt */
    inline void inc_t() { t += dt; }

    /** getter for t */
    inline double get_t() const { return t; }

    /** setter for t */
    inline void set_t(double _t) { t=_t; }

    /** to perform some second order corrections, an effective \f$ \alpha \f$ is computed here with a piecewise formula */
    inline double calc_alpha_eff(const double alpha,const double h) const
        {
        double a_eff = alpha;
        const double r = 0.1;//what is that constant, where does it come from ?
        const double M = 2.*alpha*r/dt;

        if (h>0.){ 
            if (h>M) a_eff = alpha+dt/2.*M;
            else a_eff = alpha+dt/2.*h;
            }
        else{
            if (h<-M) a_eff = alpha/(1.+dt/(2.*alpha)*M);
            else a_eff = alpha/(1.-dt/(2.*alpha)*h);
            }
        return a_eff;
        }
    
    private:
        double t;/**< physical current time of the simulation */

        double dt;/**< time-step */

        /** initialize variable time step with geometric average of DTMIN and DTMAX */
        inline void set_dt_init() {set_dt( sqrt(DTMIN*DTMAX) );}

        /** prefactor setter */
        inline void set_prefactor() {prefactor = (1.+ dt/TAUR*abs(log(dt/TAUR)));}
    };



#endif
