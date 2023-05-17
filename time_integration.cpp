#include <cfloat>
#include <ctime>
#include <signal.h>

#include "fem.h"
#include "time_integration.h"
#include "chronometer.h"
#include "fmm_demag.h"
#include "linear_algebra.h"
#include "log-stats.h"

/** Logic for finding a reasonable time step. */
class TimeStepper
    {
    const double hard_min; /**< never go below **/
    const double hard_max; /**< never exceed **/
    double soft_max;       /**< too large for this specific time step */
public:
    /** Create a TimeStepper initialized with the given initial and maximum time steps. */
    TimeStepper(double initial, double min, double max)
        : hard_min(min), hard_max(max * (1 + FLT_EPSILON)), soft_max(initial)
        {
        }

    /** Update soft_max to be no larger than max. */
    void set_soft_limit(double max) { soft_max = std::min(soft_max, max); }

    /** Return a reasonable time step. `stride' is distance to the next
     * time we want to hit exactly. */
    double operator()(double stride)
        {
        double step = std::min(stride, soft_max);
        if (step > stride - 2 * hard_min && step < stride)
            {
            step = stride - 2 * hard_min;
            }
        soft_max = std::max(soft_max, std::min(step * 1.1, hard_max));
        return step;
        }
    };

/** Summary statistics on the time steps. */
struct Stats
    {
    LogStats good_dt;    /**< dt of successful steps */
    LogStats good_dumax; /**< dumax of successful steps */
    LogStats bad_dt;     /**< dt of failed steps */
    };

static void print_stats(const Stats &s)
    {
    puts("\nTime step statistics:\n");
    puts("    time steps       count       dt [*]          dumax [*]");
    puts("    ──────────────────────────────────────────────────────────");
    printf("    successful   %9g", (double) s.good_dt.count());
    if (s.good_dt.count() != 0)
        printf("   %8.2e ± %4.2f   %8.2e ± %4.2f\n", s.good_dt.mean(), s.good_dt.stddev(),
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

/** Periodically show the percentage of work done, together with an
 * ASCII-art spinner. End the output with CR in order to keep the cursor
 * on the same line, which then gets overwritten on the next update. */
static void show_progress(double fraction_done)
    {
    const char spinner[] = "|/-\\";
    const std::chrono::duration<double> min_period(0.5);
    static int spinner_pos = 0;
    static std::chrono::time_point<std::chrono::steady_clock> last_update;
    static bool done = false;

    if (done) return;
    if (fraction_done == 1)
        {
        // Erase spinner and move to the next line.
        puts("progress:          100.00%  ");
        done = true;
        return;
        }

    auto now = std::chrono::steady_clock::now();
    if (now - last_update < min_period) return;
    last_update = now;

    printf("progress:          %.2f%% %c\r", 100 * fraction_done, spinner[spinner_pos]);
    spinner_pos = (spinner_pos + 1) % 4;
    fflush(stdout);
    }

/** Check whether we received a signal politely asking us to terminate.
 * If so, then save the current state and exit. */
static void exit_if_signal_received(const Fem &fem, const Settings &settings, const timing &t_prm,
                                    const Stats &stats)
    {
    extern volatile sig_atomic_t received_signal;  // set by signal_handler() in main.cpp
    if (!received_signal) return;
    const char *signal_name = received_signal == SIGINT    ? "SIGINT"
                              : received_signal == SIGTERM ? "SIGTERM"
                                                           : "signal";
    std::cout << "\nReceived " << signal_name;

    if (settings.save_period > 0)
        {
        std::cout << ": saving the magnetization configuration...\n";
        std::string fileName =
                settings.r_path_output_dir + '/' + settings.getSimName() + "_at_exit.sol";
        std::string metadata = settings.solMetadata(t_prm.get_t(), "idx\tmx\tmy\tmz\tphi");
        fem.msh.savesol(settings.getPrecision(), fileName, metadata);
        std::cout << "Magnetization configuration saved to " << fileName << "\n";
        }
    else
        {
        std::cout << ": magnetization configuration not saved.\n";
        }
    std::cout << "Terminating.\n";

    print_stats(stats);
    exit(1);
    }

/** compute all quantitites at time t */
inline void compute_all(Fem &fem, Settings &settings, scal_fmm::fmm &myFMM, double t)
    {
    myFMM.calc_demag(fem.msh, settings);
    fem.energy(t, settings);
    fem.evolution();
    }

int time_integration(Fem &fem, Settings &settings /**< [in] */, LinAlgebra &linAlg /**< [in] */,
                     scal_fmm::fmm &myFMM /**< [in] */, timing &t_prm, int &nt)
    {
    compute_all(fem, settings, myFMM, t_prm.get_t());

    std::string baseName = settings.r_path_output_dir + '/' + settings.getSimName();
    std::string str = baseName + ".evol";

    std::ofstream fout(str);

    if (fout.fail())
        {
        std::cout << "cannot open file " << str << std::endl;
        SYSTEM_ERROR;
        }

    fout << settings.evolMetadata(date());

    int flag(0);
    int nt_output(0);  // visible iteration count
    int status(0);     // exit status
    double t_step = settings.time_step;
    TimeStepper stepper(t_prm.get_dt(), t_prm.DTMIN, t_prm.DTMAX);
    Stats stats;

    // Loop over the visible time steps, i.e. those that will appear on the output file.
    nt = 0;
    for (double t_target = t_prm.get_t(); t_target < t_prm.tf + t_step / 2; t_target += t_step)
        {
        // Loop over the integration time steps within a visible step.
        while (t_prm.get_t() < t_target)
            {
            exit_if_signal_received(fem, settings, t_prm, stats);

            t_prm.set_dt(stepper(t_target - t_prm.get_t()));
            bool last_step = (t_prm.get_dt() == t_target - t_prm.get_t());

            if (settings.verbose)
                {
                std::cout << std::string(64, '-') << '\n';  // separator
                if (flag)
                    std::cout << "  TRYING AGAIN with a smaller time step: "
                              << "retry " << flag << '\n';
                std::cout << "evol step = " << nt_output << ", step = " << nt
                          << ", t = " << t_prm.get_t() << ", dt = " << t_prm.get_dt() << '\n';
                }

            if (t_prm.is_dt_TooSmall())
                {
                std::cout << "\n**ABORTED**: dt < DTMIN\n";
                status = 1;
                goto bailout;
                }

            Pt::pt3D Hext = settings.getValue(t_prm.get_t());

            linAlg.prepareElements(Hext, t_prm);
            int err = linAlg.solver(t_prm, nt);
            fem.vmax = linAlg.get_v_max();

            if (err)
                {
                flag++;
                stepper.set_soft_limit(t_prm.get_dt() / 2);
                stats.bad_dt.add(t_prm.get_dt());
                continue;
                }

            double dumax = t_prm.get_dt() * fem.vmax;
            stats.good_dt.add(t_prm.get_dt());
            stats.good_dumax.add(dumax);
            if (settings.verbose)
                {
                std::cout << "  -> dumax = " << dumax << ",  vmax = " << fem.vmax << std::endl;
                }

            stepper.set_soft_limit(settings.DUMAX / fem.vmax / 2);
            if (dumax > settings.DUMAX)
                {
                flag++;
                continue;
                }

            compute_all(fem, settings, myFMM, t_prm.get_t());
            nt++;
            flag = 0;

            // Prevent rounding errors from making us miss the target.
            if (last_step)
                t_prm.set_t(t_target);
            else
                t_prm.inc_t();

            if (settings.recenter) fem.recenter(settings.threshold, settings.recentering_direction);
            if (!settings.verbose) show_progress(t_prm.get_t() / t_prm.tf);
            }  // endwhile
        fem.saver(settings, t_prm, fout, nt_output++);
        }                                       // end for
    if (!settings.verbose) show_progress(1.0);  // show we are done

bailout:
    fout.close();
    print_stats(stats);
    return status;
    }
