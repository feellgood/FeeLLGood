#include <algorithm>
#include <ctime>

#include "config.h"
#include "fem.h"
#include "mesh.h"
#include "time_integration.h"

#include "tags.h"
#include "chronometer.h"

using namespace std;

void Fem::saver(Settings &settings, timing const &t_prm, ofstream &fout, const int nt) const
    {
    int save_period = settings.save_period;

    for (unsigned int i = 0; i < settings.evol_columns.size(); i++)
        {
        std::string sep;
        const std::string &keyVal = settings.evol_columns[i];
        if (i == settings.evol_columns.size() - 1)
            {
            sep = "\n";
            }
        else
            {
            sep = "\t";
            }

        if (keyVal == "iter")
            {
            fout << nt << sep;
            }
        if (keyVal == "t")
            {
            fout << t_prm.get_t() << sep;
            }
        if (keyVal == "dt")
            {
            fout << t_prm.get_dt() << sep;
            }
        if (keyVal == "max_dm")
            {
            fout << vmax * t_prm.get_dt() << sep;
            }
        if (keyVal == "<Mx>")
            {
            fout << msh.avg(Nodes::get_u_comp, Pt::IDX_X) << sep;
            }
        if (keyVal == "<My>")
            {
            fout << msh.avg(Nodes::get_u_comp, Pt::IDX_Y) << sep;
            }
        if (keyVal == "<Mz>")
            {
            fout << msh.avg(Nodes::get_u_comp, Pt::IDX_Z) << sep;
            }
        if (keyVal == "<dMx/dt>")
            {
            fout << msh.avg(Nodes::get_v_comp, Pt::IDX_X) << sep;
            }
        if (keyVal == "<dMy/dt>")
            {
            fout << msh.avg(Nodes::get_v_comp, Pt::IDX_Y) << sep;
            }
        if (keyVal == "<dMz/dt>")
            {
            fout << msh.avg(Nodes::get_v_comp, Pt::IDX_Z) << sep;
            }
        if (keyVal == "E_ex")
            {
            fout << E_exch << sep;
            }
        if (keyVal == "E_aniso")
            {
            fout << E_aniso << sep;
            }
        if (keyVal == "E_demag")
            {
            fout << E_demag << sep;
            }
        if (keyVal == "E_zeeman")
            {
            fout << E_zeeman << sep;
            }
        if (keyVal == "E_tot")
            {
            fout << Etot << sep;
            }
        if (keyVal == "Hx")
            {
            fout << settings.getField(t_prm.get_t()).x() << sep;
            }
        if (keyVal == "Hy")
            {
            fout << settings.getField(t_prm.get_t()).y() << sep;
            }
        if (keyVal == "Hz")
            {
            fout << settings.getField(t_prm.get_t()).z() << sep;
            }
        }
    fout << std::flush;

    string baseName = settings.r_path_output_dir + '/' + settings.getSimName();

    if (save_period && (nt % save_period) == 0)
        {
        string str = baseName + "_iter" + to_string(nt) + ".sol";

        if (settings.verbose)
            {
            cout << " " << str << endl;
            }

        string metadata = settings.solMetadata(t_prm.get_t(), "idx\tmx\tmy\tmz\tphi");
        msh.savesol(settings.getPrecision(), str, metadata);
        if (settings.verbose)
            {
            cout << "all nodes written." << endl;
            }
        }
    }

void Mesh::mesh::savesol(const int precision, const std::string fileName,
                         std::string const &metadata) const
    {
    ofstream fout(fileName, ios::out);
    if (fout.fail())
        {
        std::cout << "cannot open file " << fileName << std::endl;
        SYSTEM_ERROR;
        }

    fout << tags::sol::rw_time << ' ' <<  date() << '\n' << metadata << std::scientific
         << std::setprecision(precision);

    for (unsigned int i = 0; i < node.size(); i++)
        {
        fout << i << '\t' << node[i].u << '\t' << node[i].phi << endl;
        }

    fout.close();
    }

bool Mesh::mesh::savesol(const int precision, const std::string fileName,
                         std::string const &metadata, std::vector<double> const &val) const
    {
    ofstream fout(fileName, ios::out);
    if (fout.fail())
        {
        std::cout << "cannot open file " << fileName << std::endl;
        SYSTEM_ERROR;
        }

    fout << tags::sol::rw_time << ' ' << date() << '\n' << metadata << std::scientific 
         << std::setprecision(precision);

    if (node.size() == val.size())
        {
        for (unsigned int i = 0; i < node.size(); i++)
            {
            fout << i << '\t' << val[i] << endl;
            }
        }
    else
        {
        std::cout << "error: size mismatch while saving " << fileName << std::endl;
        exit(1);
        }
    fout.close();
    return !(fout.good());
    }
