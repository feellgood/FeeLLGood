#include <cfloat>
#include <signal.h>
#include "fem.h"
#include "time_integration.h"

#include "linear_algebra.h"
#include "fmm_demag.h"
#include "log-stats.h"

/** Logic for finding a reasonable time step. */
class TimeStepper {
    const double hard_min;  /**< never go below **/
    const double hard_max;  /**< never exceed **/
    double soft_max;        /**< too large for this specific time step */
public:
    /** Create a TimeStepper initialized with the given initial and maximum time steps. */
    TimeStepper(double initial, double min, double max)
    : hard_min(min), hard_max(max * (1 + FLT_EPSILON)), soft_max(initial) {}

    /** Update soft_max to be no larger than max. */
    void set_soft_limit(double max) {
        soft_max = std::min(soft_max, max);
    }

    /** Return a reasonable time step. `stride' is distance to the next
     * time we want to hit exactly. */
    double operator()(double stride) {
        double step = std::min(stride, soft_max);
        if (step > stride - 2*hard_min && step < stride) {
            step = stride - 2*hard_min;
        }
        soft_max = std::max(soft_max, std::min(step * 1.1, hard_max));
        return step;
    }
};

/** Summary statistics on the time steps. */
struct Stats {
    LogStats good_dt;  /**< dt of successful steps */
    LogStats good_dumax;  /**< dumax of successful steps */
    LogStats bad_dt;  /**< dt of failed steps */
};

static void print_stats(const Stats &s)
{
puts("\nTime step statistics:\n");
puts("    time steps       count       dt [*]          dumax [*]");
puts("    ──────────────────────────────────────────────────────────");
printf("    successful   %9g", (double) s.good_dt.count());
if (s.good_dt.count() != 0)
    printf("   %8.2e ± %4.2f   %8.2e ± %4.2f\n",
        s.good_dt.mean(), s.good_dt.stddev(),
        s.good_dumax.mean(), s.good_dumax.stddev());
else
    putchar('\n');
printf("    failed       %9g", (double) s.bad_dt.count());
if (s.bad_dt.count() != 0)
    printf("   %8.2e ± %4.2f\n\n", s.bad_dt.mean(), s.bad_dt.stddev());
else
    puts("\n");
puts("    [*] ranges given as (geometric mean) ± (relative stddev)");
}

/** Check whether we received a signal politely asking us to terminate.
 * If so, then save the current state and exit. */
static void exit_if_signal_received(const Fem &fem, const Settings &settings,
        const timing &t_prm, const Stats &stats)
{
extern volatile sig_atomic_t received_signal;  // set by signal_handler() in main.cpp
if (!received_signal) return;
const char *signal_name = received_signal==SIGINT ? "SIGINT" :
        received_signal==SIGTERM ? "SIGTERM" : "signal";
std::cout << "Received " << signal_name
        << ": saving the magnetization configuration...\n";
std::string fileName = settings.r_path_output_dir + '/'
    + settings.getSimName() + "_at_exit.sol";
fem.msh.savesol(fileName, t_prm, settings.getScale());
std::cout << "Magnetization configuration saved to "
    << fileName << "\nTerminating.\n";
print_stats(stats);
exit(1);
}

/** compute all quantitites at time t */
inline void compute_all(Fem &fem,Settings &settings,scal_fmm::fmm &myFMM,double t)
    {
    myFMM.calc_demag(fem.msh,settings);
    fem.energy(t,settings);
    fem.evolution();
    }

int time_integration(Fem &fem,Settings &settings /**< [in] */,LinAlgebra &linAlg /**< [in] */,scal_fmm::fmm &myFMM  /**< [in] */,timing &t_prm)
{
fem.DW_z  = 0.0;

compute_all(fem,settings,myFMM,t_prm.get_t());

std::string baseName = settings.r_path_output_dir + '/' + settings.getSimName();
std::string str = baseName + ".evol";

std::ofstream fout(str);

if (fout.fail())
    { std::cout << "cannot open file "<< str << std::endl; SYSTEM_ERROR; }

if(settings.evol_header)
    {
    fout << "# ";
    for(unsigned int i=0; i < (settings.evol_columns.size()-1);i++) {fout << settings.evol_columns[i] << "\t";}
    fout << settings.evol_columns[settings.evol_columns.size()-1] << "\n" << std::flush;
    }

int flag  = 0;
int nt_output = 0;  // visible iteration count
int nt = 0;         // total iteration count
double t_step = settings.time_step;
TimeStepper stepper(t_prm.get_dt(), t_prm.DTMIN, t_prm.DTMAX);
Stats stats;

// Loop over the visible time steps,
// i.e. those that will appear on the output file.
for (double t_target = t_prm.get_t(); t_target <  t_prm.tf+t_step/2; t_target += t_step)
    {
    // Loop over the integration time steps within a visible step.
    while (t_prm.get_t() < t_target)
        {
        exit_if_signal_received(fem, settings, t_prm, stats);

        t_prm.set_dt( stepper(t_target - t_prm.get_t()) );
        bool last_step = (t_prm.get_dt() == t_target - t_prm.get_t());

        if(settings.verbose)
            {
            std::cout << std::string(64, '-') << '\n';  // separator
            if (flag)
                std::cout << "  TRYING AGAIN with a smaller time step: "
                          << "retry " << flag << '\n';
            std::cout << "evol step = " << nt_output
                      << ", step = " << nt
                      << ", t = " << t_prm.get_t()
                      << ", dt = " << t_prm.get_dt() << '\n';
            }
        if (t_prm.is_dt_TooSmall()) { fem.msh.reset();break; }

        /* changement de referentiel */
        fem.DW_vz += fem.DW_dir*fem.msh.avg(Nodes::get_v_comp,Pt::IDX_Z)*fem.msh.l.z()/2.;
        
        linAlg.set_DW_vz(fem.DW_vz);
        
        Pt::pt3D Hext = settings.getValue(t_prm.get_t());
        
        linAlg.prepareElements(Hext,t_prm);
        int err = linAlg.solver(t_prm,nt);  
        fem.vmax = linAlg.get_v_max();
        
        if (err)
            {
            flag++;
            stepper.set_soft_limit(t_prm.get_dt() / 2);
            stats.bad_dt.add(t_prm.get_dt());
            continue;
            }

        double dumax = t_prm.get_dt()*fem.vmax;
        stats.good_dt.add(t_prm.get_dt());
        stats.good_dumax.add(dumax);
        if(settings.verbose) { std::cout << "  -> dumax = " << dumax << ",  vmax = "<< fem.vmax << std::endl; }
                
        stepper.set_soft_limit(settings.DUMAX / fem.vmax / 2);
        if (dumax > settings.DUMAX)
            { flag++; continue;}

        compute_all(fem,settings,myFMM,t_prm.get_t());
        
        nt++; flag=0;

        // Prevent rounding errors from making us miss the target.
        if (last_step)
            t_prm.set_t(t_target);
        else
            t_prm.inc_t();

        fem.DW_vz0 = fem.DW_vz;/* mise a jour de la vitesse du dernier referentiel et deplacement de paroi */ 
        fem.DW_z  += fem.DW_vz*t_prm.get_dt();
        if(settings.recenter) fem.recenter(settings.threshold,settings.recentering_direction);
        }//endwhile
    fem.saver(settings,t_prm,fout,nt_output++);
    }// end for

if (t_prm.is_dt_TooSmall()) { std::cout << "\n**ABORTED**: dt < DTMIN\n"; }
        
fout.close();
print_stats(stats);
return nt;
}
