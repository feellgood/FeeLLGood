#include "fem.h"
#include "config.h"

//#include <boost/format.hpp>


using namespace std;


void Fem::saver(Settings &settings, ofstream &fout, int nt)
{
int n1 = settings.n1;
int n2 = settings.n2;

if ((nt%n1)==0) {
   // fout << boost::format("%8d %+20.10e %+20.10e %+20.10e") % nt % t % dt % dumax;
   // fout << boost::format("%+20.10e %+20.10e %+20.10e") % u_moy(fem,0) % u_moy(fem,1) % u_moy(fem,2);
   // fout << boost::format("%+20.10e %+20.10e %+20.10e %+20.10e %+20.10e %+20.10e %+20.10e")% Ee % Ea % Ed % Ez % Etot % DW_z % DW_vz<< endl;
fout << nt <<"\t" << t <<"\t" << settings.dt <<"\t" << vmax*settings.dt <<"\t";
fout << moy<U>(0) <<"\t" << moy<U>(1) <<"\t" << moy<U>(2) <<"\t";

fout << E[0] <<"\t" << E[1] <<"\t" << E[2] <<"\t" << E[3] <<"\t" << Etot <<"\t" << DW_z <<"\t" << DW_vz <<"\t" << endl;
}

if ((nt%n2)==0) 
    {
    if (settings.withVtk) savecfg_vtk(settings.getSimName(),nt,nullptr);
    savesol(settings.getSimName(),settings.getScale(),nt,nullptr);

//    saveH(fem,nt);
    }
}
