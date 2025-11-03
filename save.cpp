#include <algorithm>
#include <ctime>

#include "config.h"
#include "fem.h"
#include "mesh.h"
#include "time_integration.h"

#include "tags.h"
#include "chronometer.h"

using namespace std;
using namespace Nodes;

void Fem::saver(Settings &settings, timing const &t_prm, ofstream &fout, const int nt, std::vector<Eigen::Vector3d> &s) const
    {
    int save_period = settings.save_period;
    for (unsigned int i = 0; i < settings.evol_columns.size(); i++)
        {
        const std::string &col_name = settings.evol_columns[i];

        // Split the column name as region_name:keyVal.
        int region = -1;
        std::string::size_type colon_pos = col_name.rfind(':');
        if (colon_pos != std::string::npos)
            {
            std::string region_name = col_name.substr(0, colon_pos);
            region = settings.findTetraRegionIdx(region_name);
            if (region == -1)
                {
                std::cerr << "Error: no region named '" << region_name << "'\n";
                SYSTEM_ERROR
                }
            }
        const std::string &keyVal = (region == -1) ? col_name : col_name.substr(colon_pos + 1);

        // Separator to add after this data field.
        std::string sep;
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
        else if (keyVal == "t")
            {
            fout << t_prm.get_t() << sep;
            }
        else if (keyVal == "dt")
            {
            fout << t_prm.get_dt() << sep;
            }
        else if (keyVal == "max_dm")
            {
            fout << vmax * t_prm.get_dt() << sep;
            }
        else if (keyVal == "max_angle")
            {
            fout << msh.max_angle() << sep;
            }
        else if (keyVal == "<Mx>")
            {
            fout << msh.avg(Nodes::get_u_comp, IDX_X, region) << sep;
            }
        else if (keyVal == "<My>")
            {
            fout << msh.avg(Nodes::get_u_comp, IDX_Y, region) << sep;
            }
        else if (keyVal == "<Mz>")
            {
            fout << msh.avg(Nodes::get_u_comp, IDX_Z, region) << sep;
            }
        else if (keyVal == "<dMx/dt>")
            {
            fout << msh.avg(Nodes::get_v_comp, IDX_X, region) << sep;
            }
        else if (keyVal == "<dMy/dt>")
            {
            fout << msh.avg(Nodes::get_v_comp, IDX_Y, region) << sep;
            }
        else if (keyVal == "<dMz/dt>")
            {
            fout << msh.avg(Nodes::get_v_comp, IDX_Z, region) << sep;
            }
        else if (keyVal == "E_ex")
            {
            fout << E[EXCHANGE] << sep;
            }
        else if (keyVal == "E_aniso")
            {
            fout << E[ANISOTROPY] << sep;
            }
        else if (keyVal == "E_demag")
            {
            fout << E[DEMAG] << sep;
            }
        else if (keyVal == "E_zeeman")
            {
            fout << E[ZEEMAN] << sep;
            }
        else if (keyVal == "E_tot")
            {
            fout << Etot << sep;
            }
        else if (keyVal == "Hx")
            {
            fout << settings.getField(t_prm.get_t()).x() << sep;
            }
        else if (keyVal == "Hy")
            {
            fout << settings.getField(t_prm.get_t()).y() << sep;
            }
        else if (keyVal == "Hz")
            {
            fout << settings.getField(t_prm.get_t()).z() << sep;
            }
        else
            {
            std::cerr << "Error: invalid column name '" << keyVal << "'\n";
            SYSTEM_ERROR
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

        string metadata = settings.solMetadata(t_prm.get_t());
        msh.savesol(settings.getPrecision(), str, metadata, settings.spin_acc, s);
        if (settings.verbose)
            {
            cout << "all nodes written." << endl;
            }
        }
    }

void Mesh::mesh::savesol(const int precision, const std::string fileName,
                         std::string const &metadata, bool withSpinAcc, std::vector<Eigen::Vector3d> &s) const
    {
    ofstream fout(fileName, ios::out);
    if (fout.fail())
        {
        std::cout << "cannot open file " << fileName << std::endl;
        SYSTEM_ERROR;
        }

    fout << metadata << std::scientific << std::setprecision(precision);
    Eigen::IOFormat outputSolFmt(precision, Eigen::DontAlignCols, "\t", "\t", "", "", "", "");
    for (unsigned int i = 0; i < node.size(); i++)
        {
        const int j = node_index[i];
        const dataNode &dn = node[j].d[Nodes::NEXT];
        fout << i << '\t' << dn.u.format(outputSolFmt) << '\t' << dn.phi;
        if (withSpinAcc)
            { fout << '\t' << s[j].format(outputSolFmt); }
        fout << endl;
        }
    fout.close();
    }

