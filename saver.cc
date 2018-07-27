#include "fem.h"
#include "tiny.h"

void saver(Fem &fem, Settings &settings, ofstream &fout, int nt)
{
pair <string,int> p;
p=make_pair("save_energies",-1);   int n1 = int(settings.param[p]);
p=make_pair("take_photo",-1);      int n2 = int(settings.param[p]);

double t  = fem.t;
double dt = settings.dt;
double dumax = dt * fem.vmax;
double Ee, Ea, Ed, Ez, Etot, DW_z, DW_vz;
Ee = fem.E[0]; Ea = fem.E[1]; Ed = fem.E[2]; Ez = fem.E[3]; 
Etot = fem.Etot;
DW_z  = fem.DW_z;
DW_vz = fem.DW_vz;

if ((nt%n1)==0) {
   // fout << boost::format("%8d %+20.10e %+20.10e %+20.10e") % nt % t % dt % dumax;
   // fout << boost::format("%+20.10e %+20.10e %+20.10e") % u_moy(fem,0) % u_moy(fem,1) % u_moy(fem,2);
   // fout << boost::format("%+20.10e %+20.10e %+20.10e %+20.10e %+20.10e %+20.10e %+20.10e")% Ee % Ea % Ed % Ez % Etot % DW_z % DW_vz<< endl;
fout << nt <<"\t" << t <<"\t" << dt <<"\t" << dumax <<"\t";
fout << moy<U>(fem,0) <<"\t" << moy<U>(fem,1) <<"\t" << moy<U>(fem,2) <<"\t";
fout << Ee <<"\t" << Ea <<"\t" << Ed <<"\t" << Ez <<"\t" << Etot <<"\t" << DW_z <<"\t" << DW_vz <<"\t" << endl;

}

if ((nt%n2)==0) {
//    savecfg_gnuplot(fem,nt);
//    savecfg_tecplot(fem,nt);
    if (settings.withVtk) savecfg_vtk(fem,settings.getSimName(),settings.getScale(),nt,nullptr);
    savesol(fem,settings.getSimName(),settings.getScale(),nt,nullptr);

#ifdef STAT
    savestat(fem,nt);
#endif
//    saveH(fem,nt);
    }
}
