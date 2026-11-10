#include <algorithm>
#include <ctime>

#include "config.h"
#include "fem.h"
#include "mesh.h"
#include "time_integration.h"

#include "chronometer.h"

using namespace Nodes;

void Fem::saver(Settings &settings, const timing &t_prm, std::ofstream &fout, const int nt,
        const spinAcc &spinAcc_solver) const
    {
    int save_period = settings.save_period;

    // Cache the applied field if we are going to print any of "Hx", "Hy" or "Hz".
    bool field_needed = std::any_of(settings.evol_columns.begin(), settings.evol_columns.end(),
        [](const std::string &col_name)
        { return col_name[0] == 'H' && col_name[1] >= 'x' && col_name[1] <= 'z'; });
    Eigen::Vector3d applied_field(NAN, NAN, NAN);
    if (field_needed)
        { applied_field = settings.getField(t_prm.get_t_phys()); }

    for (unsigned int i = 0; i < settings.evol_columns.size(); i++)
        {
        const std::string &col_name = settings.evol_columns[i];

        // Split the column name as region_name:keyVal.
        int region = -1;
        std::string::size_type colon_pos = col_name.rfind(':');
        if (colon_pos != std::string::npos)
            {
            std::string region_name = col_name.substr(0, colon_pos);
            region = settings.findRegionIdx<Tetra::prm>(region_name);
            if (region == -1)
                {
                std::cerr << "Error: no region named '" << region_name << "'\n";
                SYSTEM_ERROR
                }
            }
        const std::string &keyVal = (region == -1) ? col_name : col_name.substr(colon_pos + 1);

        if (keyVal == "iter")
            {
            fout << std::to_string(nt);  // not affected by std::showpos
            }
        else if (keyVal == "t")
            {
            fout << t_prm.get_t_phys();
            }
        else if (keyVal == "dt")
            {
            fout << t_prm.get_dt_phys();
            }
        else if (keyVal == "max_dm")
            {
            fout << vmax * t_prm.get_dt();
            }
        else if (keyVal == "max_angle")
            {
            fout << msh.max_angle();
            }
        else if (keyVal == "<Mx>")
            {
            fout << msh.avg(Nodes::get_u_comp, IDX_X, region);
            }
        else if (keyVal == "<My>")
            {
            fout << msh.avg(Nodes::get_u_comp, IDX_Y, region);
            }
        else if (keyVal == "<Mz>")
            {
            fout << msh.avg(Nodes::get_u_comp, IDX_Z, region);
            }
        else if (keyVal == "<dMx/dt>")
            {
            fout << gamma0 * msh.avg(Nodes::get_v_comp, IDX_X, region);
            }
        else if (keyVal == "<dMy/dt>")
            {
            fout << gamma0 * msh.avg(Nodes::get_v_comp, IDX_Y, region);
            }
        else if (keyVal == "<dMz/dt>")
            {
            fout << gamma0 * msh.avg(Nodes::get_v_comp, IDX_Z, region);
            }
        else if (keyVal == "Thiele_length")
            {
            fout << msh.thiele(region);
            }
        else if (keyVal == "E_ex")
            {
            fout << E[EXCHANGE];
            }
        else if (keyVal == "E_DMI_bulk")
            {
            fout << E[DMI_BULK];
            }
        else if (keyVal == "E_aniso")
            {
            fout << E[ANISOTROPY];
            }
        else if (keyVal == "E_demag")
            {
            fout << E[DEMAG];
            }
        else if (keyVal == "E_zeeman")
            {
            fout << E[ZEEMAN];
            }
        else if (keyVal == "E_tot")
            {
            fout << Etot;
            }
        else if (keyVal == "Hx")
            {
            fout << applied_field.x();
            }
        else if (keyVal == "Hy")
            {
            fout << applied_field.y();
            }
        else if (keyVal == "Hz")
            {
            fout << applied_field.z();
            }
        else if (keyVal == "jn_multiplier")
            { fout << settings.jnMultiplier(t_prm.get_t_phys()); }
        else
            {
            std::cerr << "Error: invalid column name '" << keyVal << "'\n";
            SYSTEM_ERROR
            }

        // Separator to add after this data field.
        fout << ((i == settings.evol_columns.size() - 1) ? '\n' : '\t');
        }
    fout << std::flush;

    std::string baseName = settings.r_path_output_dir + '/' + settings.getSimName();

    if (save_period && (nt % save_period) == 0)
        {
        std::string str = baseName + "_iter" + std::to_string(nt) + ".sol";

        if (settings.verbose)
            {
            std::cout << " " << str << "\n";
            }

        std::string metadata = settings.solMetadata(t_prm.get_t_phys());
        // potential at the time of the file (it follows the current density)
        std::vector<double> V;
        if (settings.spin_acc)
            { V = spinAcc_solver.potential(settings.jnMultiplier(t_prm.get_t_phys())); }
        msh.savesol(settings.getPrecision(), str, metadata, settings.spin_acc, spinAcc_solver.s, V);
        if (settings.verbose)
            {
            std::cout << "all nodes written.\n";
            }
        }
    }

void Mesh::mesh::savesol(const int precision, const std::string& fileName,
        const std::string &metadata, bool withSpinAcc, const std::vector<Eigen::Vector3d> &s,
        const std::vector<double> &V) const
    {
    std::ofstream fout(fileName, std::ios::out);
    if (fout.fail())
        {
        std::cerr << "cannot open file " << fileName << "\n";
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
            { fout << '\t' << s[j].format(outputSolFmt) << '\t' << V[j]; }
        fout << "\n";
        }
    fout.close();
    }
