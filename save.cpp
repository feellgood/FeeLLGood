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
        else if (keyVal == "<Mx>")
            {
            fout << msh.avg(Nodes::get_u_comp, IDX_X) << sep;
            }
        else if (keyVal == "<My>")
            {
            fout << msh.avg(Nodes::get_u_comp, IDX_Y) << sep;
            }
        else if (keyVal == "<Mz>")
            {
            fout << msh.avg(Nodes::get_u_comp, IDX_Z) << sep;
            }
        else if (keyVal == "<dMx/dt>")
            {
            fout << msh.avg(Nodes::get_v_comp, IDX_X) << sep;
            }
        else if (keyVal == "<dMy/dt>")
            {
            fout << msh.avg(Nodes::get_v_comp, IDX_Y) << sep;
            }
        else if (keyVal == "<dMz/dt>")
            {
            fout << msh.avg(Nodes::get_v_comp, IDX_Z) << sep;
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

    fout << metadata << std::scientific << std::setprecision(precision);

    Eigen::IOFormat outputSolFmt(precision, Eigen::DontAlignCols, "\t", "\t", "", "", "", "");
    for (unsigned int i = 0; i < node.size(); i++)
        {
        const dataNode &dn = node[node_index[i]].d[Nodes::NEXT];
        fout << i << '\t' << dn.u.format(outputSolFmt) << '\t' << dn.phi << endl;
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

    fout << metadata << std::scientific << std::setprecision(precision);

    if (node.size() == val.size())
        {
        for (unsigned int i = 0; i < node.size(); i++)
            {
            fout << i << '\t' << val[node_index[i]] << endl;
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
