#include <iostream>
#include <unistd.h> // for getpid()

#include "Utils/FTic.hpp"//for counter from scalfmm

#include "linear_algebra.h"
#include "fmm_demag.h"

using namespace std;

void prompt(void)
{
std::cout << "\n\t ┌────────────────────────────┐\n";
std::cout <<   "\t │         FeeLLGood          │\n";
std::cout <<   "\t │        version " << feellgood_version << "       │\n";
std::cout <<   "\t │      cnrs Grenoble-INP     │\n";
std::cout <<   "\t └────────────────────────────┘\n";
std::cout << "\t process\t\t" << getpid() << std::endl;
}



int main(int argc,char* argv[])
{
Settings mySettings;

FTic counter;
 
string fileJson;

if(argc<2)
	{
	fileJson = "settings.json";
	cout << "using default settings from " << fileJson << " JSON file." <<endl; 
	}
else 
	{
	fileJson = argv[1]; // argv[0] is "./feellgood"
	cout << "using loaded settings from " << fileJson << " JSON file." <<endl;	
	}

prompt();
	
mySettings.read(fileJson);
mySettings.infos();

Fem fem = Fem(mySettings);

counter.tic();

fem.infos();

//once fem containers are ok, linAlgebra object is built
LinAlgebra linAlg(mySettings,fem.node,fem.tet,fem.fac);
linAlg.set_Hext(fem.Hext[0],fem.Hext[1],fem.Hext[2]);

scal_fmm::fmm myFMM = scal_fmm::fmm(fem,mySettings.verbose,mySettings.scalfmmNbTh);


if(mySettings.verbose) { std::cout << "\n>> ScalFMM initialized, using " << mySettings.scalfmmNbTh << " threads.\n" << std::endl; }  

double dt0= mySettings.dt;

string baseName = mySettings.r_path_output_dir + mySettings.getSimName();
string str = baseName + ".evol";

ofstream fout(str);
if (!fout) { cerr << "cannot open file "<< str << endl; SYSTEM_ERROR; }

myFMM.calc_demag(fem,mySettings);
   
fem.DW_z  = 0.0;
fem.energy(mySettings); 
fem.evolution();

int flag  = 0;
double dt = dt0;
double t= fem.t;

int nt = 0;

fem.saver(mySettings,fout,nt);

while (t < mySettings.tf)
    {
    cout << "\n ------------------------------\n";
    if (flag) cout << "    t  : same (" << flag << ")" << endl;
    else cout << "nt = " << nt << ", t = " << t << endl; // *ct*
    cout << "dt = " << dt << endl << endl;
    if (dt < mySettings.DTMIN) { fem.reset();break; }

    /* changement de referentiel */
    fem.DW_vz += fem.DW_dir*fem.avg(Nodes::get_v_comp,Pt::IDX_Z)*fem.l.z()/2.;
    
    
    linAlg.set_DW_vz(fem.DW_vz);
    linAlg.set_dt(dt);
    int err = linAlg.solver(nt);  
    fem.vmax = linAlg.get_v_max();
    
    if (err)
        { cout << "err : " << err << endl;flag++; dt*= 0.5; mySettings.dt=dt; continue;}

    double dumax = dt*fem.vmax;
    cout << "\t dumax = " << dumax << ",  vmax = "<< fem.vmax << endl;
    if (dumax < mySettings.DUMIN) break; 
            
    if (dumax > mySettings.DUMAX)
        { flag++; dt*= 0.5; mySettings.dt=dt; continue;}
       
    myFMM.calc_demag(fem,mySettings);
       
    fem.energy(mySettings);
    if (fem.evol > 0.0)
        { cout << "Warning energy increases! : " << fem.evol << endl; }

    fem.DW_vz0 = fem.DW_vz;/* mise a jour de la vitesse du dernier referentiel et deplacement de paroi */ 
    fem.DW_z  += fem.DW_vz*dt;
    fem.evolution(); t+=dt; fem.t=t; nt++; flag=0;
    
    if(mySettings.recenter)
        {
        switch(mySettings.recentering_direction)
            {
            case 'X':fem.recentrage( mySettings.threshold,Pt::IDX_X);break;
            case 'Y':fem.recentrage( mySettings.threshold,Pt::IDX_Y);break;
            case 'Z':fem.recentrage( mySettings.threshold,Pt::IDX_Z);break;
            default: std::cout << "unknown recentering direction"<< std::endl; break;
            }
        
        }
    fem.saver(mySettings,fout,nt);
    dt = min(1.1*dt, mySettings.DTMAX); 
    mySettings.dt=dt;
    }//endwhile

if (dt < mySettings.DTMIN)
    { cout << " aborted:  dt < DTMIN"; }
        
fem.saver(mySettings,fout,nt);
        
counter.tac();
cout << "\n  * iterations: " << nt << "\n  * total computing time: " << counter.elapsed() << " s\n--- the end ---\n" << endl;
fout.close();
    
return 0;
}
