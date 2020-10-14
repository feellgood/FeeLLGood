#include <cfloat>
#include "fem.h"
#include "linear_algebra.h"
#include "fmm_demag.h"
#include "time_integration.h"

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

void compute_all(Fem &fem,Settings &settings,scal_fmm::fmm &myFMM,double t)
{
myFMM.calc_demag(fem,settings);
fem.energy(t,settings);
fem.evolution();
}

int time_integration(Fem &fem,Settings &settings /**< [in] */,LinAlgebra &linAlg /**< [in] */,scal_fmm::fmm &myFMM  /**< [in] */,timing &t_prm)
{
fem.DW_z  = 0.0;

compute_all(fem,settings,myFMM,t_prm.t);

std::string baseName = settings.r_path_output_dir + settings.getSimName();
std::string str = baseName + ".evol";

std::ofstream fout(str);
if (!fout) { std::cerr << "cannot open file "<< str << std::endl; SYSTEM_ERROR; }

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
TimeStepper stepper(t_prm.dt, t_prm.DTMIN, t_prm.DTMAX);

// Loop over the visible time steps,
// i.e. those that will appear on the output file.
for (double t_target = t_prm.t; t_target <  t_prm.tf+t_step/2; t_target += t_step)
    {
    // Loop over the integration time steps within a visible step.
    while (t_prm.t < t_target)
        {
        t_prm.dt = stepper(t_target - t_prm.t);
        bool last_step = (t_prm.dt == t_target - t_prm.t);

        if(settings.verbose)
            {
            std::cout << " ------------------------------\n";
            if (flag) std::cout << "    t  : same (" << flag << ")";
            else std::cout << "nt_output = " << nt_output << ", nt = " << nt << ", t = " << t_prm.t;
            std::cout << ", dt = " << t_prm.dt ;
            }
        if (t_prm.dt < t_prm.DTMIN) { fem.reset();break; }

        /* changement de referentiel */
        fem.DW_vz += fem.DW_dir*fem.avg(Nodes::get_v_comp,Pt::IDX_Z)*fem.l.z()/2.;
        
        linAlg.set_DW_vz(fem.DW_vz);
        //int err = linAlg.monoThreadSolver(t_prm,nt);
        Pt::pt3D Hext = settings.getValue(t_prm.t);
        int err = linAlg.solver(Hext,t_prm,nt);  
        fem.vmax = linAlg.get_v_max();
        
        if (err)
            {
            std::cout << "solver warning #" << err << ": you might need to adapt time step, max(du), and/or the refreshing period of the preconditionner" << std::endl;
            flag++;
            stepper.set_soft_limit(t_prm.dt / 2);
            continue;
            }

        double dumax = t_prm.dt*fem.vmax;
        if(settings.verbose) { std::cout << "\t dumax = " << dumax << ",  vmax = "<< fem.vmax << std::endl; }
        if (dumax < settings.DUMIN) break; 
                
        stepper.set_soft_limit(settings.DUMAX / fem.vmax / 2);
        if (dumax > settings.DUMAX)
            { flag++; continue;}

        compute_all(fem,settings,myFMM,t_prm.t);
        
        nt++; flag=0;

        // Prevent rounding errors from making us miss the target.
        if (last_step)
            t_prm.t = t_target;
        else
            t_prm.t += t_prm.dt;

        fem.DW_vz0 = fem.DW_vz;/* mise a jour de la vitesse du dernier referentiel et deplacement de paroi */ 
        fem.DW_z  += fem.DW_vz*t_prm.dt;
        if(settings.recenter) fem.recenter(settings.threshold,settings.recentering_direction);
            
        }//endwhile
    fem.saver(settings,t_prm,fout,nt_output++);
    }// end for

if (t_prm.dt < t_prm.DTMIN) { std::cout << " aborted:  dt < DTMIN"; }
        
fout.close();
return nt;
}
