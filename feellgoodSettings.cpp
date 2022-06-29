#include <iostream>
#include <string>
#include <stdio.h>  // for perror()
#include <unistd.h>  // for sysconf(), stat()
#include <sys/stat.h>  // for mkdir(), stat()
#include <sys/types.h>  // for mkdir(), stat()

#include "feellgoodSettings.h"


/***********************************************************************
 * Private helper functions.
 */

// Create the output directory if it does not exist yet.
static void create_dir_if_needed(std::string dirname)
{
    std::cout << "Output directory: " << dirname << " ";
    const char *name = dirname.c_str();
    struct stat statbuf;
    int res = stat(name, &statbuf);
    if (res != 0 && errno != ENOENT) {
        std::cout << "could not be searched.\n";
        perror(name);
        exit(1);
    }
    if (res == 0) {  // path exists
        if (S_ISDIR(statbuf.st_mode)) {
            std::cout << "(already exists)\n";
            return;
        } else {
            std::cout << "exists and is not a directory.\n";
            exit(1);
        }
    }

    // The directory does not exist (stat() reported ENOENT), create it.
    res = mkdir(name, 0777);
    if (res != 0) {
        std::cout << "could not be created.\n";
        perror(name);
        exit(1);
    }
    std::cout << "(created)\n";
}

// Bail out on errors.
static void error(const char *message)
{
    std::cerr << "CONFIGURATION ERROR: " << message << "\n";
    exit(1);
}

// Conditionally assign a variable if the node is defined.
template <typename T>
static bool assign(T &var, const YAML::Node &node)
{
    if (node) {
        var = node.as<T>();
        return true;
    }
    return false;
}

// Overload of the previous template for a unit vector.
static bool assign(Pt::pt3D &var, const YAML::Node &node)
{
    if (node) {
        if (!node.IsSequence())
            error("vectors should be YAML sequences.");
        if (node.size() != 3)
            error("vectors should have three components.");
        var = Pt::pt3D(
                node[0].as<double>(),
                node[1].as<double>(),
                node[2].as<double>()
                );
        var.normalize();
        return true;
    }
    return false;
}


/***********************************************************************
 * Public API.
 */

void Settings::infos()
{
std::cout << "\t simulation name : "<< simName << std::endl;
std::cout << "\t mesh file name: " << pbName << " with scaling factor " << getScale() << std::endl;

std::cout << "\t save energy values every " << time_step << " seconds" << std::endl;
std::cout << "\t snapshot of the magnetization configuration every " << save_period << " time steps" << std::endl;
if (restoreFileName != "")
    { std::cout << "\t restore initial magnetization distribution from file : " << restoreFileName <<std::endl; }
else
    { std::cout << "\t initial magnetization distribution from math expression" << std::endl; }

std::cout << "\t applied field Bext = [ " << sBx << ",\t" << sBy << ",\t" << sBz << " ] A/m" << std::endl;

if (recenter)
    { std::cout << "recentering along " << recentering_direction << " with threshold = " << threshold << std::endl; }

for(unsigned int i=0;i<paramTetra.size();i++) {paramTetra[i].infos();}
for(unsigned int i=0;i<paramFacette.size();i++) {paramFacette[i].infos();}

std::cout << solverNbTh+1 << " threads for assembling matrix." << std::endl;
}

void Settings::read(YAML::Node yaml)
{
    YAML::Node outputs = yaml["outputs"];
    if (outputs) {
        if (!outputs.IsMap())
            error("outputs should be a map.");
        if (assign(r_path_output_dir, outputs["directory"])) {
            // Normalize directory name.
            if (r_path_output_dir.empty())
                { r_path_output_dir = "."; }
            if (r_path_output_dir.length() > 1
                    && r_path_output_dir.back() == '/')
                { r_path_output_dir.pop_back(); }
            create_dir_if_needed(r_path_output_dir);
        }
        assign(simName, outputs["file_basename"]);
        assign(withVtk, outputs["vtk_file"]);
        assign(time_step, outputs["evol_time_step"]);
        assign(save_period, outputs["take_photo"]);
        YAML::Node columns = outputs["evol_columns"];
        if (columns) {
            if (!columns.IsSequence())
                error("outputs.evol_columns should be a sequence.");
            evol_columns.clear();
            for (auto it = columns.begin(); it != columns.end(); ++it)
                evol_columns.push_back(it->as<std::string>());
        }
        assign(evol_header, outputs["evol_header"]);
    }  // outputs

    YAML::Node mesh = yaml["mesh"];
    if (mesh) {
        if (!mesh.IsMap())
            error("mesh should be a map.");
        assign(pbName, mesh["filename"]);
        if (assign(_scale, mesh["scaling_factor"]) && _scale <= 0)
            error("mesh.scaling_factor should be positive.");
        YAML::Node volumes = mesh["volume_regions"];
        if (volumes) {
            if (!volumes.IsMap())
                error("mesh.volume_regions should be a map.");
            for (auto it = volumes.begin(); it != volumes.end(); ++it) {
                std::string name = it->first.as<std::string>();
                YAML::Node volume = it->second;
                Tetra::prm p;
                p.reg = stoi(name);
                assign(p.A, volume["Ae"]);
                assign(p.J, volume["Js"]);
                assign(p.K, volume["K"]);
                assign(p.uk, volume["uk"]);
                assign(p.K3, volume["K3"]);
                assign(p.ex, volume["ex"]);
                assign(p.ey, volume["ey"]);
                assign(p.ez, volume["ez"]);
                if (!Pt::isOrthogonal(p.ex, p.ey, p.ez, USER_TOL))
                    std::cout << "Warning: (ex, ey, ez) is not orthogonal.\n";
                assign(p.alpha_LLG, volume["alpha_LLG"]);

                // XXX: These should come from the configuration file.
                p.p_STT.beta = 0.0;
                p.p_STT.N0 = 1.0;
                p.p_STT.sigma = 1.0;
                p.p_STT.lJ = 1.0;
                p.p_STT.lsf = 1.0;
                p.p_STT.gamma0 = 1.0;
                p.p_STT.func = [](Pt::pt3D){ return 1; };

                paramTetra.push_back(p);
            }
        }  // mesh.volume_regions
        YAML::Node surfaces = mesh["surface_regions"];
        if (surfaces) {
            if (!surfaces.IsMap())
                error("mesh.surface_regions should be a map.");
            for (auto it = surfaces.begin(); it != surfaces.end(); ++it) {
                std::string name = it->first.as<std::string>();
                YAML::Node surface = it->second;
                Facette::prm p;
                p.reg = stoi(name);
                assign(p.suppress_charges, surface["suppress_charges"]);
                assign(p.Ks, surface["Ks"]);
                assign(p.uk, surface["uk"]);
                paramFacette.push_back(p);
            }
        }  // mesh.surface_regions
    }  // mesh

    YAML::Node magnetization = yaml["initial_magnetization"];
    if (magnetization) {
        if (magnetization.IsScalar()) {
            restoreFileName = magnetization.as<std::string>();
        } else if (magnetization.IsSequence()) {
            if (magnetization.size() != 3)
                error("initial_magnetization should have three components.");
            sMx = magnetization[0].as<std::string>();
            sMy = magnetization[1].as<std::string>();
            sMz = magnetization[2].as<std::string>();
            doCompile3Dprm();
        } else {
            error("initial_magnetization should be either a file name "
                  "or a vector of expressions.");
        }
    }  // initial_magnetization

    YAML::Node recentering = yaml["recentering"];
    if (recentering) {
        recenter = true;
        std::string direction = recentering["direction"].as<std::string>("Z");
        if (direction.length() != 1 || direction[0] < 'X' || direction[0] > 'Z')
            error("recentering.direction should be X, Y or Z.");
        switch (direction[0]) {
            case 'X': recentering_direction = Pt::IDX_X; break;
            case 'Y': recentering_direction = Pt::IDX_Y; break;
            case 'Z': recentering_direction = Pt::IDX_Z; break;
        }
        assign(threshold, recentering["threshold"]);
    } else {
        recenter = false;
    }  // recentering

    YAML::Node field = yaml["Bext"];
    if (field) {
        if (!field.IsSequence())
            error("Bext should be a vector of expressions.");
        if (field.size() != 3)
            error("Bext should have three components.");
        sBx = field[0].as<std::string>();
        sBy = field[1].as<std::string>();
        sBz = field[2].as<std::string>();
        doCompile1Dprm();
    }  // Bext

    YAML::Node stt = yaml["spin_transfer_torque"];
    if (stt) {
        assign(p_stt.gamma0, stt["gamma0"]);
        assign(p_stt.sigma, stt["sigma"]);
        assign(p_stt.N0, stt["dens_state"]);
        assign(p_stt.beta, stt["beta"]);
        assign(p_stt.lJ, stt["l_J"]);
        assign(p_stt.lsf, stt["l_sf"]);
        assign(p_stt.reg, stt["volume_region_reference"]);
        p_stt.func = [](Pt::pt3D){ return 1; };
        p_stt.bc = Tetra::boundary_conditions::Undef;
    } else {
        stt_flag = false;
    }  // spin_transfer_torque

    // The number of available processors (actually, hardware threads)
    // is the default for the number of threads to spin.
    int available_cpu_count = sysconf(_SC_NPROCESSORS_ONLN);

    YAML::Node solver = yaml["finite_element_solver"];
    if (solver) {
        assign(solverNbTh, solver["nb_threads"]);
        if (solverNbTh <= 0)
            solverNbTh = available_cpu_count;
        assign(MAXITER, solver["max(iter)"]);
        assign(REFRESH_PRC, solver["refresh_preconditioner_every"]);
    }  // finite_element_solver

    solver = yaml["demagnetizing_field_solver"];
    if (solver) {
        assign(scalfmmNbTh, solver["nb_threads"]);
        if (scalfmmNbTh <= 0)
            scalfmmNbTh = available_cpu_count;
    }  // demagnetizing_field_solver

    YAML::Node time_integration = yaml["time_integration"];
    if (time_integration) {
        assign(DUMAX, time_integration["max(du)"]);
        assign(dt_min, time_integration["min(dt)"]);
        assign(dt_max, time_integration["max(dt)"]);
        tf = dt_max;
        assign(tf, time_integration["final_time"]);
    }  // time_integration
}

void Settings::read(std::string filename)
{
    read(YAML::LoadFile(filename));
}
