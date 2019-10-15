#include "linear_algebra.h"
#include "fmm_demag.h"

int time_integration(Fem &fem,Settings &settings /**< [in] */,LinAlgebra &linAlg /**< [in] */,scal_fmm::fmm &myFMM  /**< [in] */)
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

double dt0= settings.dt;
int flag  = 0;
double dt = dt0;
double my_t= fem.t;
int nt_output = 0;  // visible iteration count
int nt = 0;         // total iteration count
double t_step = settings.time_step;
double reachable_dt = dt;

// Loop over the visible time steps,
// i.e. those that will appear on the output file.
for (double t_target = fem.t; t_target <  settings.tf+t_step/2; t_target += t_step)
    {
    // Loop over the integration time steps within a visible step.
    while (my_t < t_target)
        {
        std::cout << " ------------------------------\n";
        if (flag) std::cout << "    t  : same (" << flag << ")";
        else std::cout << "nt_output = " << nt_output << ", nt = " << nt << ", t = " << my_t;
        std::cout << ", dt = " << dt ;
        if (dt < settings.DTMIN) { fem.reset();break; }

        /* changement de referentiel */
        fem.DW_vz += fem.DW_dir*fem.avg(Nodes::get_v_comp,Pt::IDX_Z)*fem.l.z()/2.;
        
        linAlg.set_DW_vz(fem.DW_vz);
        linAlg.set_dt(dt);
        int err = linAlg.solver(nt);  
        fem.vmax = linAlg.get_v_max();
        
        if (err)
            { std::cout << "err : " << err << std::endl;flag++; dt*= 0.5; settings.dt=dt; continue;}

        double dumax = dt*fem.vmax;
        std::cout << "\t dumax = " << dumax << ",  vmax = "<< fem.vmax << std::endl;
        if (dumax < settings.DUMIN) break; 
                
        if (dumax > settings.DUMAX)
            { flag++; dt*= 0.5; settings.dt=dt; continue;}

        // Adjust the final steps to exactly reach the target time.
        bool last_step = false;
        reachable_dt = dt;
        if (my_t + dt >= t_target) {
            dt = t_target - my_t;
            settings.dt = dt;
            last_step = true;
        } else if (my_t + dt + settings.DTMIN > t_target) {
            dt = (t_target - my_t) / 2;
            settings.dt = dt;
        }

        myFMM.calc_demag(fem,settings);
           
        fem.energy(settings);
        if (fem.evol > 0.0)
            { std::cout << "Warning energy increases! : " << fem.evol << std::endl; }

        fem.DW_vz0 = fem.DW_vz;/* mise a jour de la vitesse du dernier referentiel et deplacement de paroi */ 
        fem.DW_z  += fem.DW_vz*dt;
        fem.evolution(); nt++; flag=0;

        // Prevent rounding errors from making us miss the target.
        if (last_step)
            my_t = t_target;
        else
            my_t += dt;
        fem.t = my_t;

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
        dt = std::min(1.1*dt, settings.DTMAX); 
        settings.dt=dt;
        }//endwhile
    if (dt < 1.1 * reachable_dt)
        {
        dt = std::min(1.1 * reachable_dt, settings.DTMAX);
        settings.dt = dt;
        }
    fem.saver(settings,fout,nt_output++);
    }// end for

if (dt < settings.DTMIN) { std::cout << " aborted:  dt < DTMIN"; }
        
fout.close();
return nt;
}
