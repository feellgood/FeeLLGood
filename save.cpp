#include <algorithm>

#include "fem.h"
#include "mesh.h"
#include "config.h"
#include "time_integration.h"

//#include <boost/format.hpp>
using namespace std;

void Fem::saver(Settings & settings,timing const& t_prm, ofstream &fout,const int nt) const
{
int save_period = settings.save_period;

// fout << boost::format("%8d %+20.10e %+20.10e %+20.10e") % nt % t % dt % dumax;
// fout << boost::format("%+20.10e %+20.10e %+20.10e") % u_moy(fem,0) % u_moy(fem,1) % u_moy(fem,2);
// fout << boost::format("%+20.10e %+20.10e %+20.10e %+20.10e %+20.10e %+20.10e %+20.10e")% Ee % Ea % Ed % Ez % Etot % DW_z % DW_vz<< endl;

for(unsigned int i = 0;i<settings.evol_columns.size();i++)
    {std::string sep;
    const std::string & keyVal = settings.evol_columns[i];
    if(i == settings.evol_columns.size() - 1) {sep = "\n";} else {sep = "\t";}

    if(keyVal == "iter") { fout << nt << sep;}
    if(keyVal == "t") { fout << t_prm.get_t() << sep;}
    if(keyVal == "dt") { fout << t_prm.get_dt() << sep;}
    if(keyVal == "max_dm") { fout << vmax*t_prm.get_dt() << sep;}
    if(keyVal == "<Mx>") { fout << msh.avg(Nodes::get_u_comp,Pt::IDX_X) << sep;}
    if(keyVal == "<My>") { fout << msh.avg(Nodes::get_u_comp,Pt::IDX_Y) << sep;}
    if(keyVal == "<Mz>") { fout << msh.avg(Nodes::get_u_comp,Pt::IDX_Z) << sep;}
    if(keyVal == "<dMx/dt>") { fout << msh.avg(Nodes::get_v_comp,Pt::IDX_X) << sep;}
    if(keyVal == "<dMy/dt>") { fout << msh.avg(Nodes::get_v_comp,Pt::IDX_Y) << sep;}
    if(keyVal == "<dMz/dt>") { fout << msh.avg(Nodes::get_v_comp,Pt::IDX_Z) << sep;}
    if(keyVal == "E_ex") { fout << E_exch << sep;}
    if(keyVal == "E_aniso") { fout << E_aniso << sep;}
    if(keyVal == "E_demag") { fout << E_demag << sep;}
    if(keyVal == "E_zeeman") { fout << E_zeeman << sep;}
    if(keyVal == "E_tot") { fout << Etot << sep;}
    if(keyVal == "DW_z") { fout << DW_z << sep;}
    if(keyVal == "DW_dz") { fout << DW_vz <<  sep;}
    if(keyVal == "Hx") { fout << settings.getValue(t_prm.get_t()).x() << sep;}
    if(keyVal == "Hy") { fout << settings.getValue(t_prm.get_t()).y() << sep;}
    if(keyVal == "Hz") { fout << settings.getValue(t_prm.get_t()).z() << sep;}
    }
fout << std::flush;

string baseName = settings.r_path_output_dir + '/' + settings.getSimName();
    
if (save_period && (nt%save_period)==0)
    {
    if (settings.withVtk)
        {
        string str = baseName + "_iter" + to_string(nt) + ".vtk";
        msh.savecfg_vtk(settings,t_prm,str);
        }
    
    string str = baseName + "_iter" + to_string(nt) + ".sol";
 
    if(settings.verbose) { cout << " " << str << endl; }    
    msh.savesol(str,t_prm,settings.getScale());
    if(settings.verbose) { cout << "all nodes written." << endl; }
    //    saveH(str);
    }
}

void Mesh::mesh::savecfg_vtk(Settings const& settings,timing const& t_prm,const string fileName) const
{
    const int TET = tet.size();
    const int NOD = node.size();

if(settings.verbose) { cout <<"\n -------------------\n " << fileName << endl; }

ofstream fout(fileName, ios::out);
if (fout.fail())
    {
    std::cout << "cannot open file : " << fileName << std::endl;
    SYSTEM_ERROR;
    }

fout << "# vtk DataFile Version 2.0" << endl;
fout << "time : " << t_prm.get_t() << endl; // boost::format(" time : %+20.10e") % fem.t;
fout << "ASCII" << endl;
fout << "DATASET UNSTRUCTURED_GRID" << endl;
fout << "POINTS "<< NOD << " float" << endl;

for(int i=0;i<NOD;i++) {fout << node[i].p << endl;}

fout << "CELLS " << setw(8) << TET << "\t" << setw(8) << (5*TET) << endl;

std::for_each(tet.begin(),tet.end(),[&fout](Tetra::Tet const &te)
    {
    fout << setw(8) << Tetra::N;    
    for (int i=0; i<Tetra::N; i++) { fout << setw(8) << te.ind[i]; }
    fout << endl;    
    });

fout << "CELL_TYPES " << setw(8) << TET << endl;
for (int t=0; t<TET; t++) { fout << setw(8) << 10 << endl; }

fout <<"POINT_DATA " << setw(8) << NOD << endl;
fout <<"SCALARS phi float 1" << endl;
fout << "LOOKUP_TABLE my_table" << endl;

for(int i=0;i<NOD;i++) {fout << node[i].phi << endl;}

fout << "VECTORS u float" << endl;
for(int i=0;i<NOD;i++) {fout << node[i].u << endl;}
//fout << boost::format("%+20.10e %+20.10e %+20.10e") % u1 % u2 % u3 << endl;
}

void Mesh::mesh::savesol(const string fileName,timing const& t_prm, const double s) const
{
ofstream fout(fileName, ios::out);
if (fout.fail())
    {
    std::cout << "cannot open file " << fileName << std::endl;
    SYSTEM_ERROR;
    }
//fout << boost::format("#time : %+20.10e ") % fem.t << endl;
fout << "#time : " << t_prm.get_t() << endl;

// fout << boost::format("%8d %+20.10f %+20.10f %+20.10f %+20.10f %+20.10f %+20.10f %+20.10e") % i % x % y % z % u1 % u2 % u3 % phi << endl;}

for(unsigned int i=0;i<node.size();i++)
        { 
        Pt::pt3D p = node[i].p/s;
        fout << i << "\t" << p << "\t" << node[i].u << "\t" << node[i].phi << '\t' << node[i].V << endl;
        }

fout.close();
}

void Mesh::mesh::saveH(const string fileName,const double t,const double scale) const
{
std::cout << " " << fileName <<"\n -------------------\n" << std::endl;

ofstream fout(fileName, ios::out);
if (fout.fail())
    {
    std::cout << "cannot open file " << fileName << std::endl;
    SYSTEM_ERROR;
    }
fout << "#time : "<< t << endl;

int idx_tet=0;
std::for_each(tet.begin(),tet.end(),[&idx_tet,&fout,scale](Tetra::Tet const &te) 
    {
    
    /*------------------- INTERPOL --------------------*/
    Pt::pt3D gauss[Tetra::NPI];
    double Hdx[Tetra::NPI], Hdy[Tetra::NPI], Hdz[Tetra::NPI];

    te.interpolation(Nodes::get_p,gauss);
    te.interpolation(Nodes::get_phi,Hdx,Hdy,Hdz);
    /*---------------------------------------------------*/
    
    for (int npi=0; npi<Tetra::NPI; npi++)
        {
	    Pt::pt3D p = gauss[npi]/scale;
	    
        fout << idx_tet << " " << npi << " " << p << " "<< Hdx[npi] << " " << Hdy[npi] << " " << Hdz[npi] << endl;
        // " %8d %3d %+20.10f %+20.10f %+20.10f %+20.10e %+20.10e %+20.10e"
        }
    idx_tet++;
    }
);//end for_each

fout.close();
}
