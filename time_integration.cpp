#include "fem.h"
#include "linear_algebra.h"
#include "fmm_demag.h"
#include "time_integration.h"

/** Logic for finding a reasonable time step. */
class TimeStepper {
    const double hard_max;  /**< never exceed **/
    double soft_max;        /**< too large for this specific time step */
public:
    /** Create a TimeStepper initialized with the given initial and maximum time steps. */
    TimeStepper(double initial, double max) : hard_max(max), soft_max(2*initial) {}

    /** Update soft_max to be no larger than max. */
    void set_soft_limit(double max) {
        soft_max = std::min(soft_max, max);
    }

    /** Return a reasonable time step. `stride' is distance to the next
     * time we want to hit exactly. */
    double operator()(double stride) {
        double target_dt = soft_max/2;
        int steps = ceil(stride / target_dt);
        soft_max = std::min(soft_max * 1.1, hard_max);
        return stride / steps;
    }
};

int time_integration(Fem &fem,Settings &settings /**< [in] */,LinAlgebra &linAlg /**< [in] */,scal_fmm::fmm &myFMM  /**< [in] */,timing &t_prm)
{
fem.DW_z  = 0.0;
fem.energy(settings); 
fem.evolution();

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
double dt = t_prm.dt;
double my_t= t_prm.t;
int nt_output = 0;  // visible iteration count
int nt = 0;         // total iteration count
double t_step = settings.time_step;
TimeStepper stepper(dt, t_prm.DTMAX);

// Loop over the visible time steps,
// i.e. those that will appear on the output file.
for (double t_target = t_prm.t; t_target <  t_prm.tf+t_step/2; t_target += t_step)
    {
    // Loop over the integration time steps within a visible step.
    while (my_t < t_target)
        {
        t_prm.dt = dt = stepper(t_target - t_prm.t);
        bool last_step = (dt == t_target - t_prm.t);

        if(settings.verbose)
            {
            std::cout << " ------------------------------\n";
            if (flag) std::cout << "    t  : same (" << flag << ")";
            else std::cout << "nt_output = " << nt_output << ", nt = " << nt << ", t = " << my_t;
            std::cout << ", dt = " << dt ;
            }
        if (dt < t_prm.DTMIN) { fem.reset();break; }

        /* changement de referentiel */
        fem.DW_vz += fem.DW_dir*fem.avg(Nodes::get_v_comp,Pt::IDX_Z)*fem.l.z()/2.;
        
        linAlg.set_DW_vz(fem.DW_vz);
        //int err = linAlg.monoThreadSolver(t_prm,nt);
        int err = linAlg.solver(t_prm,nt);  
        fem.vmax = linAlg.get_v_max();
        
        if (err)
            {
            std::cout << "err : " << err << std::endl;
            flag++;
            stepper.set_soft_limit(dt);
            continue;
            }

        double dumax = dt*fem.vmax;
        if(settings.verbose) { std::cout << "\t dumax = " << dumax << ",  vmax = "<< fem.vmax << std::endl; }
        if (dumax < settings.DUMIN) break; 
                
        stepper.set_soft_limit(settings.DUMAX / fem.vmax);
        if (dumax > settings.DUMAX)
            { flag++; continue;}

        myFMM.calc_demag(fem,settings);
           
        fem.energy(settings);
        if (settings.verbose && (fem.evol > 0.0))
            { std::cout << "Warning energy increases! : " << fem.evol << std::endl; }

        fem.DW_vz0 = fem.DW_vz;/* mise a jour de la vitesse du dernier referentiel et deplacement de paroi */ 
        fem.DW_z  += fem.DW_vz*dt;
        fem.evolution(); nt++; flag=0;

        // Prevent rounding errors from making us miss the target.
        if (last_step)
            my_t = t_target;
        else
            my_t += dt;
        t_prm.t = my_t;

        if(settings.recenter)
            {
            switch(settings.recentering_direction)
                {
                case 'X':fem.recentrage( settings.threshold,Pt::IDX_X);break;
                case 'Y':fem.recentrage( settings.threshold,Pt::IDX_Y);break;
                case 'Z':fem.recentrage( settings.threshold,Pt::IDX_Z);break;
                default: std::cout << "unknown recentering direction"<< std::endl; break;
                }
            
            }
        }//endwhile
    fem.saver(settings,t_prm,fout,nt_output++);
    }// end for

if (dt < t_prm.DTMIN) { std::cout << " aborted:  dt < DTMIN"; }
        
fout.close();
return nt;
}
