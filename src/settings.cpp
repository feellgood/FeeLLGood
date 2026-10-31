#include <iostream>
#include <sstream>
#include <string>
#include <algorithm>
#include <regex>
#include <filesystem>
#include <cmath>
#include <unistd.h>  // sysconf(), gethostname(), pipe(), fork(), read(), write(), close()
#include <sys/wait.h>

#include "tags.h"
#include "chronometer.h"
#include "settings.h"

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

// Read an axis: X, Y or Z (default Z if the node is undefined). The key is used in the error
// message.
static Nodes::index readAxis(const YAML::Node &node, const std::string &key)
    {
    std::string axis = node.as<std::string>("Z");
    if (axis.length() != 1 || axis[0] < 'X' || axis[0] > 'Z')
        { error((key + " should be X, Y or Z.").c_str()); }
    switch (axis[0])
        {
        case 'X': return Nodes::IDX_X;
        case 'Y': return Nodes::IDX_Y;
        default: return Nodes::IDX_Z;
        }
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

// Read a required number of the map `node`, the key path `where` is used in the error message.
static double readRequired(YAML::Node node, const std::string &key, const std::string &where)
    {
    if (!node[key] || !node[key].IsScalar())
        { error((where + "." + key + " is required.").c_str()); }
    return node[key].as<double>();
    }

// Read the polarizer stack of a surface region (key `slonczewski`), in physical units, see
// slonczewskiInterface.h: the F1/N interface (required), the bulk of F1 (optional, folded into F1/N
// by spinAcc) and the N/F2 interface (optional). An interface is given by its resistance-area
// product AR (Ohm m^2), its asymmetry gamma and its mixing conductance Gr, Gi (S/m^2); N/F2 may
// also give the model of its transverse spin current, transverse: sink (default) or series.
// The nodes are taken by value: operator[] of a const YAML::Node throws on a missing key.
static Slonczewski::StackSettings readStack(YAML::Node node, const std::string &where)
    {
    if (!node.IsMap())
        { error((where + " should be a map.").c_str()); }
    Slonczewski::StackSettings st;
    if (!node["p"] || node["p"].IsNull())
        { error((where + ".p is required.").c_str()); }
    st.p = { node["p"][0].as<double>(), node["p"][1].as<double>(), node["p"][2].as<double>() };
    if (st.p.norm() == 0)
        { error((where + ".p should not be zero.").c_str()); }
    st.p.normalize();

    auto readInterface = [&where](YAML::Node itf, const std::string &key)
        {
        const std::string w = where + "." + key;
        if (!itf.IsMap())
            { error((w + " should be a map.").c_str()); }
        Slonczewski::Interface i;
        i.AR = readRequired(itf, "AR", w);
        i.gamma = readRequired(itf, "gamma", w);
        i.Gr = readRequired(itf, "Gr", w);
        assign(i.Gi, itf["Gi"]);
        if (i.AR <= 0 || std::abs(i.gamma) > 1 || i.Gr < 0)
            { error((w + ": AR > 0, |gamma| <= 1 and Gr >= 0 are required.").c_str()); }
        return i;
        };

    if (!node["F1N"])
        { error((where + ".F1N is required.").c_str()); }
    st.F1N = readInterface(node["F1N"], "F1N");

    YAML::Node F1 = node["F1"];
    if (F1 && !F1.IsNull())
        {
        const std::string w = where + ".F1";
        if (!F1.IsMap())
            { error((w + " should be a map.").c_str()); }
        st.hasF1 = true;
        st.sigma1 = readRequired(F1, "sigma", w);
        st.P1bulk = readRequired(F1, "P", w);
        st.lsf1 = readRequired(F1, "l_sf", w);
        assign(st.t1, F1["thickness"]);
        if (st.sigma1 <= 0 || st.lsf1 <= 0 || st.t1 <= 0 || st.P1bulk < 0 || st.P1bulk >= 1)
            { error((w + ": sigma, l_sf, thickness > 0 and 0 <= P < 1 are required.").c_str()); }
        }

    YAML::Node NF2 = node["NF2"];
    if (NF2 && !NF2.IsNull())
        {
        st.hasNF2 = true;
        st.NF2 = readInterface(NF2, "NF2");
        YAML::Node transverse = NF2["transverse"];
        if (transverse && !transverse.IsNull())
            {
            const std::string model = transverse.as<std::string>();
            if (model == "sink")
                { st.transverseSink = true; }
            else if (model == "series")
                { st.transverseSink = false; }
            else
                { error((where + ".NF2.transverse should be sink or series.").c_str()); }
            }
        }
    return st;
    }

// Read the tunnel barrier of a surface region (key `tunnel`), see tunnelJunction.h: the
// magnetization p of F1, the resistance-area product RA (Ohm m^2) in the parallel configuration,
// the TMR, the potential V of F1 (required) and the effective tunnel polarization P1 of F1
// (optional, sqrt(TMR/(2 + TMR)) by default). The nodes are taken by value, as in readStack().
static Tunnel::Settings readTunnel(YAML::Node node, const std::string &where)
    {
    if (!node.IsMap())
        { error((where + " should be a map.").c_str()); }
    Tunnel::Settings st;
    if (!node["p"] || node["p"].IsNull())
        { error((where + ".p is required.").c_str()); }
    st.p = { node["p"][0].as<double>(), node["p"][1].as<double>(), node["p"][2].as<double>() };
    if (st.p.norm() == 0)
        { error((where + ".p should not be zero.").c_str()); }
    st.p.normalize();
    st.RA = readRequired(node, "RA", where);
    st.TMR = readRequired(node, "TMR", where);
    st.V1 = readRequired(node, "V", where);
    assign(st.P1, node["P1"]);
    if (st.RA <= 0 || st.TMR < 0)
        { error((where + ": RA > 0 and TMR >= 0 are required.").c_str()); }
    // P2 = TMR/(2 + TMR)/P1 <= 1
    if (!std::isnan(st.P1) && (st.P1 > 1 || st.P1 < Tunnel::cosineCoefficient(st.TMR)))
        { error((where + ": TMR/(2 + TMR) <= P1 <= 1 is required.").c_str()); }
    return st;
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
static const char *str(const bool x) { return x ? "true" : "false"; }

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
static const std::string str(std::string s, const int level = 0)
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
    if (s.back() == '\n')
        { s.resize(s.size() - 1); }

    // Add indentation.
    replace(s, "\n", "\n" + indent);

    return s;
    }

bool isOrthogonal(const Eigen::Ref<const Eigen::Vector3d> a, const Eigen::Ref<const Eigen::Vector3d> b,
        const Eigen::Ref<const Eigen::Vector3d> c, const double precision)
    {
    bool val = (fabs(a.dot(b)) < precision);
    val &= (fabs(b.dot(c)) < precision);
    val &= (fabs(c.dot(a)) < precision);
    return val;
    }

// Compute the SHA1 checksum of a file or a string.
// If `filename` is "-", the contents of `data` is checksummed. Otherwise the named file is
// checksummed, and `data` is ignored.
// Returns: the checksum as a 40-character hex string.
// In case of failure, prints a message to stderr and returns the empty string.
static std::string sha1sum(const std::string &filename, const std::string &data = {})
    {
    // Create a pair of pipes to communicate with the sha1sum subprocess.
    // Pipes are named from the subprocess' perspective:
    //   feellgood  --[1]--input--[0]->  sha1sum  --[1]--output--[0]->  feellgood
    int input[2], output[2];
    if (pipe(input) != 0 || pipe(output) != 0)
        {
        std::cerr << "pipe: " << strerror(errno) << '\n';
        return "";
        }

    // Run the subprocess.
    pid_t pid = fork();
    if (pid < 0)
        {
        std::cerr << "fork: " << strerror(errno) << '\n';
        return "";
        }
    if (pid == 0)
        {    // We're in the child here.
        dup2(input[0], STDIN_FILENO);
        dup2(output[1], STDOUT_FILENO);
        close(input[0]);
        close(input[1]);
        close(output[1]);
        close(output[0]);
        execlp("sha1sum", "sha1sum", filename.c_str(), NULL);
        std::cerr << "execlp: " << strerror(errno) << '\n';
        exit(1);
        }

    // We're in the parent here.
    close(input[0]);
    close(output[1]);

    // If checksumming a string, feed the string data to the subprocess.
    if (filename == "-")
        for (size_t bytes_written = 0; bytes_written < data.size(); )
            {
            const char* p = data.c_str() + bytes_written;
            size_t len = data.size() - bytes_written;
            ssize_t ret = write(input[1], p, len);
            if (ret < 0)
                {
                std::cerr << "write to sha1sum: " << strerror(errno) << '\n';
                return "";
                }
            if (ret == 0)
                {
                std::cerr << "write to sha1sum: could not write any byte\n";
                return "";
                }
            bytes_written += ret;
            }
    close(input[1]);

    // Read back the subprocess output.
    char buffer[256];
    size_t bytes_read;
    for (bytes_read = 0; bytes_read < sizeof buffer; )
        {
        char* p = buffer + bytes_read;
        size_t len = sizeof buffer - bytes_read;
        ssize_t ret = read(output[0], p, len);
        if (ret < 0)
            {
            std::cerr << "read from sha1sum: " << strerror(errno) << '\n';
            return "";
            }
        if (ret == 0)
            {
            break;
            }
        bytes_read += ret;
        }
    close(output[0]);

    // Avoid zombies.
    int wstatus;
    wait(&wstatus);
    if (!WIFEXITED(wstatus) || WEXITSTATUS(wstatus) != 0)
        {
        std::cerr << "sha1sum failed.\n";
        return "";
        }

    // Sanity check.
    std::string out_text(buffer, bytes_read);
    if (!std::regex_match(out_text, std::regex("[0-9a-f]{40}  .*\n")))
        {
        std::cerr << "sha1sum: wrong output from external command: " << out_text << '\n';
        return "";
        }

    // Remove the trailing "  <filename>\n".
    return out_text.substr(0, 40);
    }

/***********************************************************************
 * Public API.
 */

Settings::Settings()
    {
    precision = 7;  // precision is 7 digits : smaller digits of node::potential phi are varying due
                    // to residual errors
    verbose = 0;
    spin_acc = false;
    field_type = UNDEF;
    read(YAML::Load(get_default_yaml()));  // load defaults
    }

void Settings::dumpDefaults() { std::cout << get_default_yaml(); }

void Settings::dumpTemplate()
    {
    std::string s = get_default_yaml();
    s = std::regex_replace(s, std::regex(" *#.*"), "");  // remove comments
    s = std::regex_replace(s, std::regex("^\n", std::regex::multiline), "");  // remove empty lines
    std::cout << s;
    }

void Settings::toYaml()
    {
    std::cout << "outputs:\n";
    std::cout << "  directory: " << r_path_output_dir << "\n";
    std::cout << "  file_basename: " << simName << "\n";
    std::cout << "  evol_time_step: " << time_step << "\n";
    std::cout << "  final_time: " << tf << "\n";
    std::cout << "  evol_columns:\n";
    for (auto it = evol_columns.begin(); it != evol_columns.end(); ++it)
        { std::cout << "    - " << *it << "\n"; }

    std::cout << "  mag_config_every: " << save_period << "\n";
    std::cout << "  metadata:\n";
    for (auto it = userMetadata.begin(); it != userMetadata.end(); ++it)
        { std::cout << "    " << it->first << ": " << str(it->second, 2) << '\n'; }

    std::cout << "mesh:\n";
    std::cout << "  filename: " << pbName << "\n";
    std::cout << "  length_unit: " << _scale << "\n";
    std::cout << "  volume_regions:\n";
    for (auto it = paramTetra.begin(); it != paramTetra.end(); ++it)
        {
        if (it->regName == "__default__")  // skip
            { continue; }
        std::cout << "    " << it->regName << ":\n";
        std::cout << "      Ae: " << it->A << "\n";
        std::cout << "      Ms: " << it->Ms << "\n";
        std::cout << "      K: " << it->K << "\n";
        if (it->K != 0)
            { std::cout << "      uk: " << str(it->uk) << "\n"; }
        std::cout << "      K3: " << it->K3 << "\n";
        if (it->K3 != 0)
            {
            std::cout << "      ex: " << str(it->e.col(IDX_X)) << "\n";
            std::cout << "      ey: " << str(it->e.col(IDX_Y)) << "\n";
            std::cout << "      ez: " << str(it->e.col(IDX_Z)) << "\n";
            }
        std::cout << "      alpha_LLG: " << it->alpha_LLG << "\n";
        std::cout << "      sigma: " << it->sigma << "\n";
        std::cout << "      dens_state: " << it->N0 << "\n";
        std::cout << "      P: " << it->P << "\n";
        std::cout << "      l_sd: " << it->lsd << "\n";
        std::cout << "      l_sf: " << it->lsf << "\n";
        }
    std::cout << "  surface_regions:\n";
    for (auto it = paramTriangle.begin(); it != paramTriangle.end(); ++it)
        {
        if (it->regName == "__default__")  // skip
            { continue; }
        std::cout << "    " << it->regName << ":\n";
        std::cout << "      suppress_charges: " << str(it->suppress_charges) << "\n";
        std::cout << "      Ks: " << it->Ks << "\n";
        if (it->Ks != 0)
            { std::cout << "      uk: " << str(it->uk) << "\n"; }
        std::cout <<  "      J:";
        if (isnan(it->jn))
            { std::cout << "\n"; }
        else if (it->jn_current_density)
            { std::cout << " current_density\n"; }
        else
            { std::cout << ' ' << it->jn << '\n'; }
        std::cout <<  "      V:";
        if (isnan(it->V))
            { std::cout << "\n"; }
        else
            { std::cout << ' ' << it->V << '\n'; }
        if (spin_acc)
            {
            std::cout << "      uP: " << str(it->uP) << "\n";
            std::cout << "      s: " << str(it->s) << "\n";
            if (it->slonczewski)
                {
                const Slonczewski::StackSettings &st = it->stack;
                auto itf = [](const Slonczewski::Interface &i)
                    {
                    std::ostringstream o;
                    o << "{AR: " << i.AR << ", gamma: " << i.gamma << ", Gr: " << i.Gr
                      << ", Gi: " << i.Gi << "}";
                    return o.str();
                    };
                std::cout << "      slonczewski:\n";
                std::cout << "        p: " << str(st.p) << "\n";
                if (st.hasF1)
                    {
                    std::cout << "        F1: {sigma: " << st.sigma1 << ", P: " << st.P1bulk
                              << ", l_sf: " << st.lsf1 << ", thickness: " << st.t1 << "}\n";
                    }
                std::cout << "        F1N: " << itf(st.F1N) << "\n";
                if (st.hasNF2)
                    {
                    std::cout << "        NF2: " << itf(st.NF2) << ", transverse: "
                              << (st.transverseSink ? "sink" : "series") << "\n";
                    }
                }
            if (it->tunnel)
                {
                const Tunnel::Settings &b = it->barrier;
                const Tunnel::Polarizations P = Tunnel::polarizations(b);
                std::cout << "      tunnel:\n";
                std::cout << "        p: " << str(b.p) << "\n";
                std::cout << "        RA: " << b.RA << "\n";
                std::cout << "        TMR: " << b.TMR << "\n";
                std::cout << "        P1: " << P.P1 << "  # P2 = " << P.P2 << "\n";
                std::cout << "        V: " << b.V1 << "\n";
                }
            }
        }
    std::cout << "initial_magnetization: ";
    if (!sM.empty())
        { std::cout << str(sM) << "\n"; }
    else if (restoreFileName.empty())
        { std::cout << "[\"" << sMx << "\", \"" << sMy << "\", \"" << sMz << "\"]\n"; }
    else
        { std::cout << restoreFileName << "\n"; }
    std::cout << "initial_time:";
    if (isnan(initial_time))
        { std::cout << '\n'; }
    else
        { std::cout << ' ' << initial_time << '\n'; }
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
    std::cout << "absorbing_ends:\n";
    std::cout << "  enable: " << str(absorb.enable) << "\n";
    if (absorb.enable)
        {
        std::cout << "  axis: " << "XYZ"[absorb.axis] << "\n";
        std::cout << "  length: " << absorb.length << "\n";
        std::cout << "  alpha_max: " << absorb.alpha_max << "\n";
        std::cout << "  regions:";
        if (absorb.regions.empty())
            { std::cout << '\n'; }
        else
            {
            std::cout << " [";
            for (unsigned int i = 0; i < absorb.regions.size(); i++)
                { std::cout << (i ? ", " : "") << absorb.regions[i]; }
            std::cout << "]\n";
            }
        }
    std::cout << "Bext: ";
    if (field_type == R4toR3)
        { std::cout << "\n  space: " << str(sB_space, 1) << "\n  time: " << str(sB_time, 1) << '\n'; }
    else if (!sB.empty())
        { std::cout << str(sB) << "\n"; }
    else
        { std::cout << "[\"" << sBx << "\", \"" << sBy << "\", \"" << sBz << "\"]\n"; }
    std::cout << "current_density:";
    if (hasCurrentDensity())
        { std::cout << ' ' << str(sCurrentDensity) << '\n'; }
    else
        { std::cout << '\n'; }
    std::cout << "spin_accumulation:\n";
    std::cout << "  enable: " << str(spin_acc) << "\n";
    if (spin_acc)
        {
        std::cout << "  anisotropic_diffusion: " << str(anisotropic_diffusion) << "\n";
        std::cout << "  V_file: " << str(V_file) << "\n";
        }
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
    ss << tags::common::version << ' ' << feellgood_version << "\n";
    char name[HOST_NAME_MAX];
    if (gethostname(name, HOST_NAME_MAX) != ENAMETOOLONG)
        { ss << tags::common::hostname << ' ' << name << "\n"; }
    ss << tags::common::rw_time << ' ' << date() << "\n";
    for (size_t i = 0; i < settings_checksums.size(); ++i)
        {
        const auto& item = settings_checksums[i];
        ss << "## settings file[" << i << "].filename: " << item.first << "\n";
        if (!item.second.empty())
            {
            ss << "## settings file[" << i << "].checksum: sha1:" << item.second << "\n";
            }
        }
    ss << "## mesh.filename: " << mesh_checksum.first << "\n";
    ss << "## mesh.checksum: sha1:" << mesh_checksum.second << "\n";

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
    ss << tags::common::columns << ' ';
    for (unsigned int i = 0; i < (evol_columns.size() - 1); i++)
        { ss << evol_columns[i] << '\t'; }
    ss << evol_columns[evol_columns.size() - 1] << "\n";
    return ss.str();
    }

std::string Settings::solMetadata(const double t) const
    {
    std::ostringstream ss = commonMetadata();
    ss << tags::sol::time << ' ' << std::scientific << t << "\n";
    ss << tags::common::columns << ' ' << tags::sol::defaultColumnsTitle;
    if (spin_acc)
        { ss << '\t' << tags::sol::sColumnsTitle; }
    ss << "\n";
    return ss.str();
    }

void Settings::read(YAML::Node yaml)
    {
    YAML::Node outputs = yaml["outputs"];
    if (outputs && !outputs.IsNull())
        {
        if (!outputs.IsMap())
            { error("outputs should be a map."); }
        if (assign(r_path_output_dir, outputs["directory"]))
            {
            // Normalize directory name.
            if (r_path_output_dir.empty())
                { r_path_output_dir = "."; }

            if (r_path_output_dir.length() > 1 && r_path_output_dir.back() == '/')
                { r_path_output_dir.pop_back(); }
            }
        assign(simName, outputs["file_basename"]);
        assign(time_step, outputs["evol_time_step"]);
        assign(tf, outputs["final_time"]);

        YAML::Node mag_config_every = outputs["mag_config_every"];
        if (mag_config_every.Scalar() == "true")  // catch an easy mistake
            { error("outputs.mag_config_every should be an integer or `false'."); }
        else if (mag_config_every.Scalar() == "false")
            { save_period = 0; }
        else if (assign(save_period, mag_config_every) && save_period < 0)
                { save_period = 0; }

        YAML::Node columns = outputs["evol_columns"];
        if (columns && !columns.IsNull())
            {
            if (!columns.IsSequence())
                { error("outputs.evol_columns should be a sequence."); }
            evol_columns.clear();
            for (auto it = columns.begin(); it != columns.end(); ++it)
                { evol_columns.push_back(it->as<std::string>()); }
            }

        YAML::Node metadata = outputs["metadata"];
        if (metadata && !metadata.IsNull())
            {
            if (!metadata.IsMap())
                { error("outputs.metadata should be a map."); }
            for (auto it = metadata.begin(); it != metadata.end(); ++it)
                {
                // Search the key within the already known user metadata.
                std::string key = it->first.as<std::string>();
                auto pos = std::find_if(userMetadata.begin(), userMetadata.end(),
                        [key](const MetadataItem &item){ return item.first == key; });
                MetadataItem item = {key, it->second.as<std::string>()};
                if (pos != userMetadata.end())  // if found, replace
                    { *pos = item; }
                else  // otherwise, append to the list
                    { userMetadata.push_back(item); }
                }
            }
        }  // outputs

    YAML::Node mesh = yaml["mesh"];
    if (mesh && !mesh.IsNull())
        {
        if (!mesh.IsMap())
            { error("mesh should be a map."); }
        assign(pbName, mesh["filename"]);
        if (assign(_scale, mesh["length_unit"]) && _scale <= 0)
            { error("mesh.length_unit should be positive."); }
        YAML::Node volumes = mesh["volume_regions"];
        if (volumes && !volumes.IsNull())
            {
            if (!volumes.IsMap())
                { error("mesh.volume_regions should be a map."); }
            int default_idx = findRegionIdx<Tetra::prm>("__default__");
            for (auto it = volumes.begin(); it != volumes.end(); ++it)
                {
                bool has_Ms = false, has_Js = false;  // whether Ms and/or Js are defined
                std::string name = it->first.as<std::string>();
                YAML::Node volume = it->second;
                Tetra::prm &p = paramTetra.emplace_back();
                if (default_idx >= 0)
                    { p = paramTetra[default_idx]; }
                p.regName = name;

                assign(p.A, volume["Ae"]);
                if (assign(p.Ms, volume["Ms"]))
                    { has_Ms = true; }
                double Js;
                if (assign(Js, volume["Js"]))
                    { has_Js = true; p.Ms = Js / mu0; }
                if (has_Ms && has_Js)
                    { error("Ms and Js should not be both defined in the same region."); }

                assign(p.K, volume["K"]);
                assign(NORMALIZE, p.uk, volume["uk"]);
                assign(p.K3, volume["K3"]);
                Eigen::Vector3d ex,ey,ez;
                assign(NORMALIZE, ex, volume["ex"]);
                assign(NORMALIZE, ey, volume["ey"]);
                assign(NORMALIZE, ez, volume["ez"]);
                p.e << ex, ey, ez;
                if (!isOrthogonal(p.e.col(IDX_X), p.e.col(IDX_Y), p.e.col(IDX_Z), USER_TOL))
                    { std::cerr << "Warning: (ex, ey, ez) is not orthogonal.\n"; }
                assign(p.alpha_LLG, volume["alpha_LLG"]);
                assign(p.sigma, volume["sigma"]);
                assign(p.N0, volume["dens_state"]);
                assign(p.P, volume["P"]);
                assign(p.lsd, volume["l_sd"]); // exists only for magnetic material, as Ms
                assign(p.lsf, volume["l_sf"]);
                }
            }  // mesh.volume_regions
        YAML::Node surfaces = mesh["surface_regions"];
        if (surfaces && !surfaces.IsNull())
            {
            if (!surfaces.IsMap())
                { error("mesh.surface_regions should be a map."); }
            int default_idx = findRegionIdx<Triangle::prm>("__default__");
            for (auto it = surfaces.begin(); it != surfaces.end(); ++it)
                {
                std::string name = it->first.as<std::string>();
                YAML::Node surface = it->second;
                Triangle::prm &p = paramTriangle.emplace_back();
                if (default_idx >= 0)
                    { p = paramTriangle[default_idx]; }
                p.regName = name;
                assign(p.suppress_charges, surface["suppress_charges"]);
                assign(p.Ks, surface["Ks"]);
                assign(NORMALIZE, p.uk, surface["uk"]);

                // J and V may be null, which we map to NAN. jn may also be the function
                // current_density(t) of the settings: jn is then the reference 1 A/m^2.
                YAML::Node jn = surface["jn"];
                p.jn_current_density = false;
                if (jn && jn.IsScalar() && jn.as<std::string>() == "current_density")
                    {
                    p.jn = 1.0;
                    p.jn_current_density = true;
                    }
                else if (!assign(p.jn, jn))
                    { p.jn = NAN; }
                if (!assign(p.V, surface["V"]))
                    { p.V = NAN; }
                if (!isnan(p.jn) && !isnan(p.V))
                    { error("A surface region cannot have both no-null J and V."); }
                if (!assign(NORMALIZE, p.uP, surface["uP"]))
                    { p.uP = Eigen::Vector3d(NAN,NAN,NAN); }
                if (!assign(!NORMALIZE, p.s, surface["s"]))
                    { p.s = Eigen::Vector3d(NAN,NAN,NAN); }
                p.slonczewski = false;
                YAML::Node stack = surface["slonczewski"];
                if (stack && !stack.IsNull())
                    {
                    const std::string where = "mesh.surface_regions." + name + ".slonczewski";
                    p.slonczewski = true;
                    p.stack = readStack(stack, where);
                    if (isnan(p.jn))
                        { error((where + " requires jn on the same surface.").c_str()); }
                    if (std::isfinite(p.s.norm()) || std::isfinite(p.uP.norm()))
                        { error((where + " excludes s and uP on the same surface.").c_str()); }
                    }
                p.tunnel = false;
                YAML::Node barrier = surface["tunnel"];
                if (barrier && !barrier.IsNull())
                    {
                    const std::string where = "mesh.surface_regions." + name + ".tunnel";
                    p.tunnel = true;
                    p.barrier = readTunnel(barrier, where);
                    if (!isnan(p.jn) || !isnan(p.V) || std::isfinite(p.s.norm())
                        || std::isfinite(p.uP.norm()) || p.slonczewski)
                        {
                        error((where + " excludes jn, V, s, uP and slonczewski on the same "
                               "surface.").c_str());
                        }
                    }
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
                { restoreFileName = s_mag; }
            else
                {
                sM = s_mag;
                mag_parser.set_function(sM);
                }
            }
        else if (magnetization.IsSequence())
            {
            if (magnetization.size() != DIM)
                { error("initial_magnetization should have three components."); }
            sMx = magnetization[0].as<std::string>();
            sMy = magnetization[1].as<std::string>();
            sMz = magnetization[2].as<std::string>();
            mag_parser.set_expressions("x,y,z", sMx, sMy, sMz);
            }
        else
            { error("initial_magnetization should be a file name or a vector of expressions."); }
        }  // initial_magnetization

    assign(initial_time, yaml["initial_time"]);

    YAML::Node recentering = yaml["recentering"];
    if (recentering && !recentering.IsNull())
        {
        assign(recenter, recentering["enable"]);
        recentering_direction = readAxis(recentering["direction"], "recentering.direction");
        assign(threshold, recentering["threshold"]);
        }  // recentering

    YAML::Node absorbing = yaml["absorbing_ends"];
    if (absorbing && !absorbing.IsNull())
        {
        if (!absorbing.IsMap())
            { error("absorbing_ends should be a map."); }
        // the section is enabled by its presence, unless enable is given
        if (!assign(absorb.enable, absorbing["enable"]))
            { absorb.enable = true; }
        absorb.axis = readAxis(absorbing["axis"], "absorbing_ends.axis");
        assign(absorb.length, absorbing["length"]);
        assign(absorb.alpha_max, absorbing["alpha_max"]);
        YAML::Node regions = absorbing["regions"];
        if (regions)
            {
            absorb.regions.clear();
            if (regions.IsScalar())  // a single region, as a sequence of one element
                { absorb.regions.push_back(regions.as<std::string>()); }
            else if (!regions.IsNull())
                {
                if (!regions.IsSequence())
                    {
                    error("absorbing_ends.regions should be a volume region or a sequence of "
                          "volume regions.");
                    }
                for (auto it = regions.begin(); it != regions.end(); ++it)
                    { absorb.regions.push_back(it->as<std::string>()); }
                }
            }
        if (absorb.enable && !(absorb.length > 0))
            { error("absorbing_ends.length should be positive."); }
        }  // absorbing_ends

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
            if (field.size() != DIM)
                { error("Bext should have three components."); }
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
                {
                error(
                "Bext should be a function or a vector of expressions or space & time expressions."
                );
                }
            }
        }  // Bext

    YAML::Node current_density = yaml["current_density"];
    if (current_density && !current_density.IsNull())
        {
        if (!current_density.IsScalar())
            { error("current_density should be a function of t or an expression of t."); }
        sCurrentDensity = current_density.as<std::string>();
        const bool is_function = std::regex_search(sCurrentDensity,
                                                   std::regex("^\\s*function\\b"));
        const std::string js_current_density = is_function
                ? sCurrentDensity : "function(t) { return (" + sCurrentDensity + "); }";
        current_density_parser.set_function(js_current_density);
        // the time functions of the applied field can call current_density(t)
        field_parser.define_global_function("current_density", js_current_density);
        field_time_parser.define_global_function("current_density", js_current_density);
        }  // current_density

    YAML::Node spAcc = yaml["spin_accumulation"];
    if (spAcc && !spAcc.IsNull())
        {
        assign(spin_acc, spAcc["enable"]);
        assign(anisotropic_diffusion, spAcc["anisotropic_diffusion"]);
        assign(V_file, spAcc["V_file"]);
	    }  // spin_accumulation

    // The number of available processors (actually, hardware threads) is the default for the number
    // of threads to spin.
    int available_cpu_count = sysconf(_SC_NPROCESSORS_ONLN);

    YAML::Node solver = yaml["demagnetizing_field_solver"];
    if (solver && !solver.IsNull())
        {
        assign(scalfmmNbTh, solver["nb_threads"]);
        if (scalfmmNbTh <= 0)
            { scalfmmNbTh = available_cpu_count; }
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
        if (pos != simName.npos)
            { simName.erase(0, pos + 1); }

        // Remove everything after last dot.
        pos = simName.rfind('.');
        if (pos != simName.npos)
            { simName.erase(pos); }
        }
    }

bool Settings::read(const std::string& filename)
    {
    std::string file_display_name;
    std::string config_src;
    YAML::Node config;
    if (filename == "-")
        {
        file_display_name = getFileDisplayName();
        std::ostringstream str_stream;
        std::cin >> str_stream.rdbuf();
        config_src = str_stream.str();
        config = YAML::Load(config_src);
        }
    else
        {
        try
            {
            config = YAML::LoadFile(filename);
            file_display_name = std::filesystem::canonical(filename);
            }
        catch (const YAML::BadFile &)
            {
            return false;
            }
        }
    settings_checksums.push_back(std::make_pair(file_display_name, sha1sum(filename, config_src)));
    if (config.IsNull())
        { return false; }
    read(config);
    if (!pbName.empty())
        {
        std::string mesh_filename = std::filesystem::canonical(pbName);
        mesh_checksum = std::make_pair(mesh_filename, sha1sum(pbName));
        }
    return true;
    }
