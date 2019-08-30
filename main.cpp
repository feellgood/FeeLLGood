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

void calc_demag(Fem &fem,Settings &mySettings,OctreeClass *tree,KernelClass *kernels)
{
fmm::demag<0, CellClass, ContainerClass, LeafClass, OctreeClass, KernelClass, FmmClass> (fem,mySettings, tree, kernels); // Hd(u)
fmm::demag<1, CellClass, ContainerClass, LeafClass, OctreeClass, KernelClass, FmmClass> (fem,mySettings, tree, kernels); // Hd(v), second order contribution
}

int main(int argc,char* argv[])
{
Settings mySettings;
Fem fem;
FTic counter;
OctreeClass *tree    = nullptr;
KernelClass *kernels = nullptr; 
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
fem.readMesh(mySettings);
counter.tic();
fem.femutil(mySettings);
fem.chapeaux(mySettings.EPSILON);

if (mySettings.restore)
    { fem.readSol(mySettings.getScale(), mySettings.restoreFileName); }
else
    {
    std::cout<< "initial magnetization(x,y,z,t=0) :\nMx =" << mySettings.sMx << "\nMy =" << mySettings.sMy << "\nMz =" << mySettings.sMz <<endl; 
    fem.init_distrib(mySettings);
    }

fem.direction(Pt::IDX_Z);/* determination de la direction de propagation de la paroi */
fem.t=0.;
fem.infos();

//once fem containers are ok, linAlgebra object is built
LinAlgebra linAlg(mySettings,fem.NOD,fem.node,fem.tet,fem.fac);

fmm::init< CellClass, ContainerClass, LeafClass, OctreeClass, KernelClass, FmmClass> (fem, mySettings.scalfmmNbTh, tree, kernels);

cout << "init scalfmm done.\n" << endl;

double dt0= mySettings.dt;

fem.Hext[0] = nu0*mySettings.Bext[0];
fem.Hext[1] = nu0*mySettings.Bext[1];
fem.Hext[2] = nu0*mySettings.Bext[2]; 

cout << "Bext : " << mySettings.Bext[0] << "\t" << mySettings.Bext[1] << "\t" << mySettings.Bext[2] << endl;

string baseName = mySettings.r_path_output_dir + mySettings.getSimName();
string str = baseName + ".evol";

ofstream fout(str);
if (!fout) { cerr << "cannot open file "<< str << endl; SYSTEM_ERROR; }

calc_demag(fem,mySettings, tree,kernels);
   
fem.DW_z  = 0.0;
fem.energy(mySettings); 
fem.evolution();

int flag  = 0;
double dt = dt0;
double t= fem.t;
fem.vmax  = 0.0;
fem.DW_vz = fem.DW_vz0 = 0.0;
fem.DW_z  = 0.0;
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
    fem.DW_vz += fem.DW_dir*fem.avg<V>(Pt::IDX_Z)*fem.l.z()/2.;
    
    linAlg.set_Hext(fem.Hext[0],fem.Hext[1],fem.Hext[2]);
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
       
    calc_demag(fem,mySettings, tree, kernels);
       
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
cout << "\n  * iterations: " << nt << "\n  * total computing time: " << counter.elapsed() << " s\n" << endl;
fout.close();
//      initialisation du temps pour l'iteration en champ suivante
fem.t=0.;
mySettings.dt=dt0;
    
cout << "--- the end ---" << endl;
delete tree;
delete kernels;

return 0;
}
