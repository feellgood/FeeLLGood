#include <algorithm>

#include "fem.h"
#include "config.h"

//#include <boost/format.hpp>
using namespace std;

void Fem::saver(Settings const& settings, ofstream &fout,const int nt) const
{
int n1 = settings.n1;
int n2 = settings.n2;

if ((nt%n1)==0)
    {
   // fout << boost::format("%8d %+20.10e %+20.10e %+20.10e") % nt % t % dt % dumax;
   // fout << boost::format("%+20.10e %+20.10e %+20.10e") % u_moy(fem,0) % u_moy(fem,1) % u_moy(fem,2);
   // fout << boost::format("%+20.10e %+20.10e %+20.10e %+20.10e %+20.10e %+20.10e %+20.10e")% Ee % Ea % Ed % Ez % Etot % DW_z % DW_vz<< endl;
    
    for(unsigned int i = 0;i<settings.evol_columns.size();i++)
        {std::string sep;
        const std::string & keyVal = settings.evol_columns[i];
        if(i == settings.evol_columns.size() - 1) {sep = "\n";} else {sep = "\t";}
    
        if(keyVal == "iter") { fout << nt << sep;}
        if(keyVal == "t") { fout << t << sep;}
        if(keyVal == "dt") { fout << settings.dt << sep;}
        if(keyVal == "max dm") { fout << vmax*settings.dt << sep;}
        if(keyVal == "<mx>") { fout << avg(Nodes::get_u_comp,Pt::IDX_X) << sep;}
        if(keyVal == "<my>") { fout << avg(Nodes::get_u_comp,Pt::IDX_Y) << sep;}
        if(keyVal == "<mz>") { fout << avg(Nodes::get_u_comp,Pt::IDX_Z) << sep;}
        if(keyVal == "E_ex") { fout << E_exch << sep;}
        if(keyVal == "E_aniso") { fout << E_aniso << sep;}
        if(keyVal == "E_demag") { fout << E_demag << sep;}
        if(keyVal == "E_zeeman") { fout << E_zeeman << sep;}
        if(keyVal == "E_tot") { fout << Etot << sep;}
        if(keyVal == "DW_z") { fout << DW_z << sep;}
        if(keyVal == "DW_dz") { fout << DW_vz <<  sep;}    
        }
    fout << std::flush;
    }

    string baseName = settings.r_path_output_dir + settings.getSimName();
    
if ((nt%n2)==0) 
    {
    if (settings.withVtk)
        {
        string str = baseName + "_iter" + to_string(nt) + ".vtk";
        savecfg_vtk(settings,str);
        }
    
    string str = baseName + "_iter" + to_string(nt) + ".sol";
 
    if(settings.verbose) { cout << " " << str << endl; }    
    savesol(str,settings.getScale());
    if(settings.verbose) { cout << "all nodes written." << endl; }
    //    saveH(str);
    }
}

void Fem::savecfg_vtk(Settings const& settings,const string fileName) const
{
    const int TET = tet.size();
    const int NOD = node.size();
    
if(settings.verbose) { cout <<"\n -------------------\n " << fileName << endl; }

ofstream fout(fileName, ios::out);
if (!fout)
    {
    if(settings.verbose) {cerr << "cannot open file : " << fileName << endl;}
    SYSTEM_ERROR;
    }

fout << "# vtk DataFile Version 2.0" << endl;
fout << "time : " << t << endl; // boost::format(" time : %+20.10e") % fem.t;
fout << "ASCII" << endl;
fout << "DATASET UNSTRUCTURED_GRID" << endl;
fout << "POINTS "<< NOD << " float" << endl;

std::for_each(node.begin(),node.end(),[this,&fout](Nodes::Node const&n)
        {fout << n.p.x() << "\t" << n.p.y() << "\t" << n.p.z() + DW_z << endl;});
//fout << boost::format("%+20.10e %+20.10e %+20.10e") % x % y % z << endl;

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

std::for_each(node.begin(),node.end(),[&fout](Nodes::Node const& n){fout << n.phi << endl;} );

fout << "VECTORS u float" << endl;
std::for_each(node.begin(),node.end(),[&fout](Nodes::Node const& n){fout << n.u << endl;} );
//fout << boost::format("%+20.10e %+20.10e %+20.10e") % u1 % u2 % u3 << endl;
}

void Fem::savesol(const string fileName, const double s) const
{
ofstream fout(fileName, ios::out);
if (!fout)
    {
    cerr << "cannot open file " << fileName << endl;
    SYSTEM_ERROR;
    }
//fout << boost::format("#time : %+20.10e ") % fem.t << endl;
fout << "#time : " << t << endl;

// fout << boost::format("%8d %+20.10f %+20.10f %+20.10f %+20.10f %+20.10f %+20.10f %+20.10e") % i % x % y % z % u1 % u2 % u3 % phi << endl;}

int i=0;
std::for_each(node.begin(),node.end(),[&i,s,&fout](Nodes::Node const&n)
        { 
        Pt::pt3D p = n.p/s;
        fout << i << "\t" << p << "\t" << n.u << "\t" << n.phi << endl;
        i++;
        } // lambda works with << overloaded
    );

fout.close();
}

void Fem::saveH(const string fileName,const double scale) const
{
cout << " " << fileName << endl <<" -------------------" << endl << endl;

ofstream fout(fileName, ios::out);
if (!fout){
   cerr << "cannot open file " << fileName << endl;
   SYSTEM_ERROR;}
fout << "#time : "<< t << endl;

int idx_tet=0;
std::for_each(tet.begin(),tet.end(),[&idx_tet,&fout,scale](Tetra::Tet const &te) 
    {
    
    /*------------------- INTERPOL --------------------*/
    double gauss[DIM][Tetra::NPI];
    double Hdx[Tetra::NPI], Hdy[Tetra::NPI], Hdz[Tetra::NPI];

    te.interpolation(Nodes::get_p,gauss);
    te.interpolation(Nodes::get_phi,Hdx,Hdy,Hdz);
    /*---------------------------------------------------*/
    
    for (int npi=0; npi<Tetra::NPI; npi++)
        {
	    double x = gauss[0][npi]/scale;
	    double y = gauss[1][npi]/scale;
	    double z = gauss[2][npi]/scale;
	
        fout << idx_tet << " " << npi << " " << x << " " << y << " " << z << " "<< Hdx[npi] << " " << Hdy[npi] << " " << Hdz[npi] << endl;
        // " %8d %3d %+20.10f %+20.10f %+20.10f %+20.10e %+20.10e %+20.10e"
        }
    idx_tet++;
    }
);//end for_each

fout.close();
}
