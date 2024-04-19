#ifndef time_integration_h
#define time_integration_h

/** all timing parameters for integrating in time LLG with adaptative time-step and relaxation corrections throught prefactor=f(dt) */
class timing
    {
public:
    /** constructor : fully initializes the timing parameters,
     initial time step initialized with geometric average of DTMIN and DTMAX */
    inline timing(const double _tf, const double _dtmin, const double _dtmax)
        : tf(_tf), DTMIN(_dtmin), DTMAX(_dtmax), TAUR(100. * DTMAX), t(0)
        {
        set_dt(sqrt(DTMIN * DTMAX));
        }

    /** final time of the simulation */
    const double tf;

    /** minimum step time for time integrator */
    const double DTMIN;  // 1e-14;

    /** maximum step time for time integrator */
    const double DTMAX;  //  1e-5 en stat ;  1e-7 en dynamique;

    /** reduced \f$ \tau_r = 100 DT_{\mathrm{max}} \f$ */
    const double TAUR;

    /** this prefactor must be synchronized with dt value */
    double prefactor;

    /** getter for time step dt */
    inline double get_dt() const { return dt; }

    /** setter for time step dt : prefactor is computed from the new dt value to maintain
     * synchronization of dt and prefactor */
    inline void set_dt(const double _dt)
        {
        dt = _dt;
        double t_tilde = _dt / TAUR;
        prefactor = (1. + t_tilde * abs(log( t_tilde)));
        }

    /** basic test to know if dt is too small */
    inline bool is_dt_TooSmall() const { return (dt < DTMIN); }

    /** increment time t with time step dt */
    inline void inc_t() { t += dt; }

    /** getter for t */
    inline double get_t() const { return t; }

    /** setter for t */
    inline void set_t(double _t) { t = _t; }

private:
    double t; /**< physical current time of the simulation */
    double dt; /**< time-step */
    };

#endif
