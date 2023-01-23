#include <algorithm>

#include "fem.h"
#include "mesh.h"
#include "config.h"
#include "time_integration.h"

using namespace std;

void Fem::saver(Settings & settings,timing const& t_prm, ofstream &fout,const int nt) const
{
int save_period = settings.save_period;

for(unsigned int i = 0;i<settings.evol_columns.size();i++)
    {
    std::string sep;
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
        msh.savecfg_vtk(t_prm,str);
        }
    
    string str = baseName + "_iter" + to_string(nt) + ".sol";
 
    if(settings.verbose)
        { cout << " " << str << endl; }

    string metadata = settings.buildMetadata(t_prm.get_t(),"idx\tx\ty\tz\tmx\tmy\tmz\tphi");
    msh.savesol(settings.getPrecision(),str,metadata,settings.getScale());
    if(settings.verbose)
        { cout << "all nodes written." << endl; }
    //saveH(str,t_prm.get_t(),settings.getScale());
    }
}

void Mesh::mesh::savecfg_vtk(timing const& t_prm,const string fileName) const
{
    const int TET = tet.size();
    const int NOD = node.size();

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

for(int i=0;i<NOD;i++)
    {fout << node[i].p << endl;}

fout << "CELLS " << setw(8) << TET << "\t" << setw(8) << (5*TET) << endl;

std::for_each(tet.begin(),tet.end(),[&fout](Tetra::Tet const &te)
    {
    fout << setw(8) << Tetra::N;    
    for (int i=0; i<Tetra::N; i++) 
        { fout << setw(8) << te.ind[i]; }
    fout << endl;    
    });

fout << "CELL_TYPES " << setw(8) << TET << endl;
for (int t=0; t<TET; t++)
    { fout << setw(8) << 10 << endl; }

fout <<"POINT_DATA " << setw(8) << NOD << endl;
fout <<"SCALARS phi float 1" << endl;
fout << "LOOKUP_TABLE my_table" << endl;

for(int i=0;i<NOD;i++)
    {fout << node[i].phi << endl;}

fout << "VECTORS u float" << endl;
for(int i=0;i<NOD;i++)
    {fout << node[i].u << endl;}
}

void Mesh::mesh::savesol(const int precision, const string fileName, string const& metadata, const double s) const
    {
    ofstream fout(fileName, ios::out);
    if (fout.fail())
        {
        std::cout << "cannot open file " << fileName << std::endl;
        SYSTEM_ERROR;
        }

    fout << metadata << std::scientific << std::setprecision(precision);

    for(unsigned int i=0;i<node.size();i++)
        { 
        Pt::pt3D p = node[i].p/s;
        fout << i << "\t" << p << "\t" << node[i].u << "\t" << node[i].phi << endl;
        }

    fout.close();
    }

bool Mesh::mesh::savesol(const int precision, const string fileName, string const& metadata, std::vector<double> const& val) const
    {
    ofstream fout(fileName, ios::out);
    if (fout.fail())
        {
        std::cout << "cannot open file " << fileName << std::endl;
        SYSTEM_ERROR;
        }

    fout << metadata << std::scientific << std::setprecision(precision);

    if (node.size() == val.size())
        {
        for(unsigned int i=0;i<node.size();i++)
            { fout << i << '\t' << val[i] << endl; }
        }
    else
        { std::cout << "error: size mismatch while saving " << fileName << std::endl; exit(1); }
    fout.close();
    
    return !(fout.good());
    }

void Mesh::mesh::saveH(const string fileName,const double t,const double scale) const
    {
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
        Pt::pt3D gauss[Tetra::NPI];
        double Hdx[Tetra::NPI], Hdy[Tetra::NPI], Hdz[Tetra::NPI];

        te.interpolation(Nodes::get_p,gauss);
        te.interpolation(Nodes::get_phi,Hdx,Hdy,Hdz);

        for (int npi=0; npi<Tetra::NPI; npi++)
            {
            Pt::pt3D p = gauss[npi]/scale;
            fout << idx_tet << " " << npi << " " << p << " "<< Hdx[npi] << " " << Hdy[npi] << " " << Hdz[npi] << endl;
            }
        idx_tet++;
        });//end for_each

    fout.close();
    }
