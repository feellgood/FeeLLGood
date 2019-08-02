#include <iostream>
#include <unistd.h> // for getpid()

#include "Utils/FTic.hpp"//for counter

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
vector<Seq> seq;
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
	
mySettings.read(fileJson,seq);
mySettings.infos(seq);
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
int nseq=0;
for (vector<Seq>::iterator it = seq.begin(); it!=seq.end(); ++it)
    {
    double &Bini=it->Bini;
    double &Bfin=it->Bfin;
    double &dB=it->dB;
    triple &a = it->a;
    nseq++;
    fem.SEQ=nseq;

    cout << "dB : " << dB << endl;

    for (int loop=0; loop<=(int)((Bfin-Bini)/dB+0.5); loop++)
    	{
    	double Bext=Bini+dB*loop;
        fem.Bext=Bext;
        fem.Hext[0]=nu0*Bext*a[0];
        fem.Hext[1]=nu0*Bext*a[1];
        fem.Hext[2]=nu0*Bext*a[2]; 

        cout << "Bext : " << Bext*a[0] << "\t" << Bext*a[1] << "\t" << Bext*a[2] << endl;

        string baseName = mySettings.r_path_output_dir + mySettings.getSimName();
        string str = baseName +"_"+ to_string(fem.SEQ) + "_B" + to_string(fem.Bext) + ".evol";

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
            fem.DW_vz += fem.DW_dir*fem.moy<V>(Pt::IDX_Z)*fem.l.z()/2.;
    
            linAlg.set_Hext(fem.Hext[0],fem.Hext[1],fem.Hext[2]);
            linAlg.set_DW_vz(fem.DW_vz);
            linAlg.set_dt(dt);
            int err = linAlg.vsolve(nt);  
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
            double mz = fem.moy<U>(Pt::IDX_Z);	    
            if(mySettings.recentering)
                { fem.recentrage( 0.1,Pt::IDX_Z,mz); }
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
        }//endfor Bext
    }//endfor it

cout << "--- the end ---" << endl;
delete tree;
delete kernels;

return 0;
}
