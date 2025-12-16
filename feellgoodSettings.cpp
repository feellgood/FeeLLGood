#include <iostream>
#include <string>
#include <algorithm>
#include <cmath>
#include <unistd.h>  // for sysconf(), gethostname()

#include "tags.h"
#include "chronometer.h"
#include "feellgoodSettings.h"

using namespace Nodes;

/***********************************************************************
 * Access to the default configuration embedded from the file
 * default-settings.yml.
 */

extern "C"
    {
    extern char _binary_default_settings_yml_start[];
    extern char _binary_default_settings_yml_end[];
    }

/***********************************************************************
 * Private helper functions.
 */

// Return the YAML document defining the defaults.
static std::string get_default_yaml()
    {
    return std::string(_binary_default_settings_yml_start,
                       _binary_default_settings_yml_end - _binary_default_settings_yml_start);
    }

// Bail out on errors.
static void error(const char *message)
    {
    std::cerr << "CONFIGURATION ERROR: " << message << "\n";
    exit(1);
    }

// Conditionally assign a variable if the node is defined and scalar.
template<typename T>
static bool assign(T &var, const YAML::Node &node)
    {
    if (node.IsScalar())
        {
        var = node.as<T>();
        return true;
        }
    return false;
    }

// Overload of the previous template for a 3D vector, with optional normalization.
static bool assign(const bool _NORMALIZE, Eigen::Vector3d &var, const YAML::Node &node)
    {
    if (node && !node.IsNull())
        {
        if (!node.IsSequence()) error("vectors should be YAML sequences.");
        if (node.size() != 3) error("vectors should have three components.");
        var = { node[0].as<double>(), node[1].as<double>(), node[2].as<double>() };
        if(_NORMALIZE) var.normalize();
        return true;
        }
    return false;
    }

// Replace, within `s', all occurrences of `search' by `replacement'.
static void replace(std::string &s, const std::string &search, const std::string &replacement)
    {
    for (size_t pos = 0; (pos = s.find(search, pos)) != s.npos;)
        {
        s.replace(pos, search.size(), replacement);
        pos += replacement.size();
        }
    }

// Stringify a boolean
static const char *str(bool x) { return x ? "true" : "false"; }

// Stringify a vector
// do NOT use std::to_string here, because cout << to_string(1e-7); does print 0.000000
// with C++20 we should use std::format
static const std::string str(Eigen::Vector3d v)
    {
    char buffer[100];
    std::sprintf(buffer,"[%g, %g, %g]",v.x(),v.y(),v.z());
    return std::string(buffer);
    }

// Stringify a string: if it contains a newline, convert it to a multiline string in "literal style"
// (introduced by '|'). Otherwise enclose it in quotes, and escape embedded quotes.
// `level` is the indentation level of the property whose value is being stringified.
static const std::string str(std::string s, int level = 0)
    {
    // If it doesn't have a newline, quote it.
    if (s.find('\n') == s.npos)
        {
        replace(s, "\"", "\\\"");  // escape embedded quotes
        return '"' + s + '"';
        }

    // Indentation to be added to each line.
    int indent_size = (level + 1) * 2;  // the value is indented one level more than the property
    std::string indent(indent_size, ' ');

    // Prepend "|\n".
    s.insert(0, "|\n");

    // Remove trailing eol.
    if (s.back() == '\n') s.resize(s.size() - 1);

    // Add indentation.
    replace(s, "\n", "\n" + indent);

    return s;
    }

bool isOrthogonal(Eigen::Ref<Eigen::Vector3d> a, Eigen::Ref<Eigen::Vector3d> b, Eigen::Ref<Eigen::Vector3d> c, const double precision)
    {
    bool val = (fabs(a.dot(b)) < precision);
    val &= (fabs(b.dot(c)) < precision);
    val &= (fabs(c.dot(a)) < precision);
    return val;
    }

/***********************************************************************
 * Public API.
 */

Settings::Settings()
    {
    precision = 7;  // precision is 7 digits : smaller digits of node::potential phi are varying due
                    // to residual errors
    verbose = 0;
    withTsv = true;
    spin_acc = false;
    field_type = UNDEF;
    read(YAML::Load(get_default_yaml()));  // load defaults
    }

void Settings::dumpDefaults() { std::cout << get_default_yaml(); }

void Settings::toYaml()
    {
    std::cout << "outputs:\n";
    std::cout << "  directory: " << r_path_output_dir << "\n";
    std::cout << "  file_basename: " << simName << "\n";
    std::cout << "  evol_time_step: " << time_step << "\n";
    std::cout << "  final_time: " << tf << "\n";
    std::cout << "  evol_columns:\n";
    for (auto it = evol_columns.begin(); it != evol_columns.end(); ++it)
        {
        std::cout << "    - " << *it << "\n";
        }
    std::cout << "  mag_config_every: " << save_period << "\n";
    std::cout << "  metadata:\n";
    for (auto it = userMetadata.begin(); it != userMetadata.end(); ++it)
        {
        std::cout << "    " << it->first << ": " << str(it->second, 2) << '\n';
        }
    std::cout << "mesh:\n";
    std::cout << "  filename: " << pbName << "\n";
    std::cout << "  length_unit: " << _scale << "\n";
    std::cout << "  volume_regions:\n";
    for (auto it = paramTetra.begin(); it != paramTetra.end(); ++it)
        {
        if (it->regName == "__default__")  // skip
            continue;
        std::cout << "    " << it->regName << ":\n";
        std::cout << "      Ae: " << it->A << "\n";
        std::cout << "      Js: " << it->J << "\n";
        std::cout << "      K: " << it->K << "\n";
        if (it->K != 0) std::cout << "      uk: " << str(it->uk) << "\n";
        std::cout << "      K3: " << it->K3 << "\n";
        if (it->K3 != 0)
            {
            std::cout << "      ex: " << str(it->ex) << "\n";
            std::cout << "      ey: " << str(it->ey) << "\n";
            std::cout << "      ez: " << str(it->ez) << "\n";
            }
        std::cout << "      alpha_LLG: " << it->alpha_LLG << "\n";
        std::cout << "      sigma: " << it->sigma << "\n";
        std::cout << "      dens_state: " << it->N0 << "\n";
        std::cout << "      P: " << it->P << "\n";
        std::cout << "      l_sd: " << it->lsd << "\n";
        std::cout << "      l_sf: " << it->lsf << "\n";
        std::cout << "      spin_hall: " << it->spinHall << "\n";
        }
    std::cout << "  surface_regions:\n";
    for (auto it = paramFacette.begin(); it != paramFacette.end(); ++it)
        {
        if (it->regName == "__default__")  // skip
            continue;
        std::cout << "    " << it->regName << ":\n";
        std::cout << "      suppress_charges: " << str(it->suppress_charges) << "\n";
        std::cout << "      Ks: " << it->Ks << "\n";
        if (it->Ks != 0) std::cout << "      uk: " << str(it->uk) << "\n";
        std::cout <<  "      J:";
        if (isnan(it->J))
            std::cout << "\n";
        else
            std::cout << ' ' << it->J << '\n';
        std::cout <<  "      V:";
        if (isnan(it->V))
            std::cout << "\n";
        else
            std::cout << ' ' << it->V << '\n';
        if (spin_acc)
            {
            std::cout << "      uP: " << str(it->uP) << "\n";
            std::cout << "      s: " << str(it->s) << "\n";
            }
        }
    std::cout << "initial_magnetization: ";
    if (!sM.empty())
        std::cout << str(sM) << "\n";
    else if (restoreFileName.empty())
        std::cout << "[\"" << sMx << "\", \"" << sMy << "\", \"" << sMz << "\"]\n";
    else
        std::cout << restoreFileName << "\n";
    std::cout << "initial_time:";
    if (isnan(initial_time))
        std::cout << '\n';
    else
        std::cout << ' ' << initial_time << '\n';
    std::cout << "recentering:\n";
    std::cout << "  enable: " << str(recenter) << "\n";
    if (recenter)
        {
        std::cout << "  direction: ";
        switch (recentering_direction)
            {
            case IDX_UNDEF: std::cout << "UNDEF\n"; break;
            case IDX_X: std::cout << "X\n"; break;
            case IDX_Y: std::cout << "Y\n"; break;
            case IDX_Z: std::cout << "Z\n"; break;
            }
        std::cout << "  threshold: " << threshold << "\n";
        }
    std::cout << "Bext: ";
    if (field_type == R4toR3)
        std::cout << "\n  space: " << str(sB_space, 1) << "\n  time: " << str(sB_time, 1) << '\n';
    else if (!sB.empty())
        std::cout << str(sB) << "\n";
    else
        std::cout << "[\"" << sBx << "\", \"" << sBy << "\", \"" << sBz << "\"]\n";
    std::cout << "spin_accumulation:\n";
    std::cout << "  enable: " << str(spin_acc) << "\n";
    if (spin_acc)
        { std::cout << "  V_file: " << str(V_file) << "\n"; }
    std::cout << "demagnetizing_field_solver:\n";
    std::cout << "  nb_threads: " << scalfmmNbTh << "\n";
    std::cout << "finite_element_solver:\n";
    std::cout << "  max(iter): " << MAXITER << "\n";
    std::cout << "  tolerance: " << TOL << "\n";
    std::cout << "time_integration:\n";
    std::cout << "  max(du): " << DUMAX << "\n";
    std::cout << "  min(dt): " << dt_min << "\n";
    std::cout << "  max(dt): " << dt_max << "\n";
    }

std::ostringstream Settings::commonMetadata() const
    {
    std::ostringstream ss;
    ss << tags::evol::version << ' ' << feellgood_version << std::endl;
    char name[HOST_NAME_MAX];
    if (gethostname(name, HOST_NAME_MAX) != ENAMETOOLONG)
        {
        ss << tags::evol::hostname << ' ' << name << std::endl;
        }
    ss << tags::evol::rw_time << ' ' << date() << std::endl;
    ss << tags::evol::settings_file << ' ' << getFileDisplayName() << std::endl;
    for (auto it = userMetadata.begin(); it != userMetadata.end(); ++it)
        {
        std::string value = it->second;
        if (value.back() == '\n') value.resize(value.size() - 1);  // remove trailing newline
        replace(value, "\n", "\n##  ");  // escape and indent continuation lines
        ss << "## " << it->first << ": " << value << '\n';
        }
    return ss;
    }

std::string Settings::evolMetadata() const
    {
    std::ostringstream ss = commonMetadata();
    ss << tags::evol::columns << ' ';
    for (unsigned int i = 0; i < (evol_columns.size() - 1); i++)
        {
        ss << evol_columns[i] << '\t';
        }
    ss << evol_columns[evol_columns.size() - 1] << std::endl;
    return ss.str();
    }

std::string Settings::solMetadata(const double t) const
    {
    std::ostringstream ss = commonMetadata();
    ss << tags::sol::time << ' ' << std::scientific << t << std::endl;
    ss << tags::sol::columns << ' ' << tags::sol::defaultColumnsTitle;
    if (spin_acc)
        ss << '\t' << tags::sol::sColumnsTitle;
    ss << std::endl;
    return ss.str();
    }

void Settings::read(YAML::Node yaml)
    {
    YAML::Node outputs = yaml["outputs"];
    if (outputs && !outputs.IsNull())
        {
        if (!outputs.IsMap()) error("outputs should be a map.");
        if (assign(r_path_output_dir, outputs["directory"]))
            {
            // Normalize directory name.
            if (r_path_output_dir.empty())
                {
                r_path_output_dir = ".";
                }
            if (r_path_output_dir.length() > 1 && r_path_output_dir.back() == '/')
                {
                r_path_output_dir.pop_back();
                }
            }
        assign(simName, outputs["file_basename"]);
        assign(time_step, outputs["evol_time_step"]);
        assign(tf, outputs["final_time"]);
        YAML::Node mag_config_every = outputs["mag_config_every"];
        if (mag_config_every.Scalar() == "true")
            {  // catch an easy mistake
            error("outputs.mag_config_every should be an integer or `false'.");
            }
        else if (mag_config_every.Scalar() == "false")
            {
            save_period = 0;
            }
        else
            {
            if (assign(save_period, mag_config_every) && save_period < 0)
                {
                save_period = 0;
                }
            }
        YAML::Node columns = outputs["evol_columns"];
        if (columns && !columns.IsNull())
            {
            if (!columns.IsSequence()) error("outputs.evol_columns should be a sequence.");
            evol_columns.clear();
            for (auto it = columns.begin(); it != columns.end(); ++it)
                evol_columns.push_back(it->as<std::string>());
            }
        YAML::Node metadata = outputs["metadata"];
        if (metadata && !metadata.IsNull())
            {
            if (!metadata.IsMap()) error("outputs.metadata should be a map.");
            for (auto it = metadata.begin(); it != metadata.end(); ++it)
                {
                // Search the key within the already known user metadata.
                std::string key = it->first.as<std::string>();
                auto pos = std::find_if(userMetadata.begin(), userMetadata.end(),
                        [key](const MetadataItem &item){ return item.first == key; });
                MetadataItem item = {key, it->second.as<std::string>()};
                if (pos != userMetadata.end())  // if found, replace
                    *pos = item;
                else  // otherwise, append to the list
                    userMetadata.push_back(item);
                }
            }
        }  // outputs

    YAML::Node mesh = yaml["mesh"];
    if (mesh && !mesh.IsNull())
        {
        if (!mesh.IsMap()) error("mesh should be a map.");
        assign(pbName, mesh["filename"]);
        if (assign(_scale, mesh["length_unit"]) && _scale <= 0)
            error("mesh.length_unit should be positive.");
        YAML::Node volumes = mesh["volume_regions"];
        if (volumes && !volumes.IsNull())
            {
            if (!volumes.IsMap()) error("mesh.volume_regions should be a map.");
            int default_idx = findTetraRegionIdx("__default__");
            for (auto it = volumes.begin(); it != volumes.end(); ++it)
                {
                std::string name = it->first.as<std::string>();
                YAML::Node volume = it->second;
                Tetra::prm &p = paramTetra.emplace_back();
                if (default_idx >= 0) p = paramTetra[default_idx];
                p.regName = name;
                assign(p.A, volume["Ae"]);
                assign(p.J, volume["Js"]);
                assign(p.K, volume["K"]);
                assign(NORMALIZE, p.uk, volume["uk"]);
                assign(p.K3, volume["K3"]);
                assign(NORMALIZE, p.ex, volume["ex"]);
                assign(NORMALIZE, p.ey, volume["ey"]);
                assign(NORMALIZE, p.ez, volume["ez"]);
                if (!isOrthogonal(p.ex, p.ey, p.ez, USER_TOL))
                    std::cout << "Warning: (ex, ey, ez) is not orthogonal.\n";
                assign(p.alpha_LLG, volume["alpha_LLG"]);

                assign(p.sigma, volume["sigma"]);
                assign(p.N0, volume["dens_state"]);
                assign(p.P, volume["P"]);
                assign(p.lsd, volume["l_sd"]); // exists only for magnetic material, as Js
                assign(p.lsf, volume["l_sf"]);
                assign(p.spinHall, volume["spin_hall"]); // SOT contribution to spin diffusion,
                                                         // through spin Hall effect
                }
            }  // mesh.volume_regions
        YAML::Node surfaces = mesh["surface_regions"];
        if (surfaces && !surfaces.IsNull())
            {
            if (!surfaces.IsMap()) error("mesh.surface_regions should be a map.");
            int default_idx = findFacetteRegionIdx("__default__");
            for (auto it = surfaces.begin(); it != surfaces.end(); ++it)
                {
                std::string name = it->first.as<std::string>();
                YAML::Node surface = it->second;
                Facette::prm &p = paramFacette.emplace_back();
                if (default_idx >= 0) p = paramFacette[default_idx];
                p.regName = name;
                assign(p.suppress_charges, surface["suppress_charges"]);
                assign(p.Ks, surface["Ks"]);
                assign(NORMALIZE, p.uk, surface["uk"]);

                // J and V may be null, which we map to NAN.
                if (!assign(p.J, surface["J"]))
                    p.J = NAN;
                if (!assign(p.V, surface["V"]))
                    p.V = NAN;
                if (!isnan(p.J) && !isnan(p.V))
                    error("A surface region cannot have both no-null J and V.");
                if (!assign(NORMALIZE, p.uP, surface["uP"]))
                    {p.uP = Eigen::Vector3d(NAN,NAN,NAN); }
                if (!assign(!NORMALIZE, p.s, surface["s"]))
                    { p.s = Eigen::Vector3d(NAN,NAN,NAN); }
                }
            }  // mesh.surface_regions
        }      // mesh

    YAML::Node magnetization = yaml["initial_magnetization"];
    if (magnetization && !magnetization.IsNull())
        {
        if (magnetization.IsScalar())
            {
            std::string s_mag = magnetization.as<std::string>();
            if (s_mag.find("function") == std::string::npos || s_mag.find('{') == std::string::npos)
                {
                restoreFileName = s_mag;
                }
            else
                {
                sM = s_mag;
                mag_parser.set_function(sM);
                }
            }
        else if (magnetization.IsSequence())
            {
            if (magnetization.size() != DIM)
                error("initial_magnetization should have three components.");
            sMx = magnetization[0].as<std::string>();
            sMy = magnetization[1].as<std::string>();
            sMz = magnetization[2].as<std::string>();
            mag_parser.set_expressions("x,y,z", sMx, sMy, sMz);
            }
        else
            {
            error("initial_magnetization should be a file name or a vector of expressions.");
            }
        }  // initial_magnetization

    assign(initial_time, yaml["initial_time"]);

    YAML::Node recentering = yaml["recentering"];
    if (recentering && !recentering.IsNull())
        {
        assign(recenter, recentering["enable"]);
        std::string direction = recentering["direction"].as<std::string>("Z");
        if (direction.length() != 1 || direction[0] < 'X' || direction[0] > 'Z')
            error("recentering.direction should be X, Y or Z.");
        switch (direction[0])
            {
            case 'X': recentering_direction = IDX_X; break;
            case 'Y': recentering_direction = IDX_Y; break;
            case 'Z': recentering_direction = IDX_Z; break;
            }
        assign(threshold, recentering["threshold"]);
        }  // recentering

    YAML::Node field = yaml["Bext"];
    if (field && !field.IsNull())
        {
        if (field.IsScalar())
            {
            sB = field.as<std::string>();
            field_parser.set_function(sB);
            field_type = RtoR3;
            }
        else if (field.IsSequence())
            {
            if (field.size() != DIM) error("Bext should have three components.");
            sBx = field[0].as<std::string>();
            sBy = field[1].as<std::string>();
            sBz = field[2].as<std::string>();
            field_parser.set_expressions("t", sBx, sBy, sBz);
            field_type = RtoR3;
            }
        else
            {
            YAML::Node field_space = field["space"];
            YAML::Node field_time = field["time"];
            if (field_space && field_time)
                {
                if (field_space.IsScalar())
                    {
                    sB_space = field_space.as<std::string>();
                    field_space_parser.set_function(sB_space);
                    }
                if (field_time.IsScalar())
                    {
                    sB_time = field_time.as<std::string>();
                    field_time_parser.set_function(sB_time);
                    }
                field_type = R4toR3;
                }
            else
                error("Bext should be a function or a vector of expressions or space & time expressions.");
            }
        }  // Bext

    YAML::Node spAcc = yaml["spin_accumulation"];
    if (spAcc && !spAcc.IsNull())
        {
        assign(spin_acc, spAcc["enable"]);
        assign(V_file, spAcc["V_file"]);
	    }  // spin_accumulation

    // The number of available processors (actually, hardware threads) is the default for the number
    // of threads to spin.
    int available_cpu_count = sysconf(_SC_NPROCESSORS_ONLN);

    YAML::Node solver = yaml["demagnetizing_field_solver"];
    if (solver && !solver.IsNull())
        {
        assign(scalfmmNbTh, solver["nb_threads"]);
        if (scalfmmNbTh <= 0) scalfmmNbTh = available_cpu_count;
        }  // demagnetizing_field_solver

    solver = yaml["finite_element_solver"];
    if (solver && !solver.IsNull())
        {
        assign(MAXITER, solver["max(iter)"]);
        assign(TOL,solver["tolerance"]);
        }  // finite_element_solver

    YAML::Node time_integration = yaml["time_integration"];
    if (time_integration && !time_integration.IsNull())
        {
        assign(DUMAX, time_integration["max(du)"]);
        assign(dt_min, time_integration["min(dt)"]);
        assign(dt_max, time_integration["max(dt)"]);
        }  // time_integration

    // outputs.file_basename defaults to base name of mesh.filename.
    if (simName.empty() && !pbName.empty())
        {
        simName = pbName;

        // Remove directory part.
        size_t pos = simName.rfind('/');
        if (pos != simName.npos) simName.erase(0, pos + 1);

        // Remove everything after last dot.
        pos = simName.rfind('.');
        if (pos != simName.npos) simName.erase(pos);
        }
    }

bool Settings::read(std::string filename)
    {
    YAML::Node config;
    if (filename == "-")
        config = YAML::Load(std::cin);
    else
        {
        try
            {
            config = YAML::LoadFile(filename);
            }
        catch (const YAML::BadFile &)
            {
            return false;
            }
        }
    if (config.IsNull()) return false;
    read(config);
    return true;
    }
