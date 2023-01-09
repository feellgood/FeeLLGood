#include <iostream>
#include <string>
#include <unistd.h>  // for sysconf()

#include "feellgoodSettings.h"


/***********************************************************************
 * Access to the default configuration embedded from the file
 * default-settings.yml.
 */

extern "C" {
    extern char _binary_default_settings_yml_start[];
    extern char _binary_default_settings_yml_end[];
}


/***********************************************************************
 * Private helper functions.
 */

// Return the YAML document defining the defaults.
static std::string get_default_yaml()
{
    return std::string(
            _binary_default_settings_yml_start,
            _binary_default_settings_yml_end
            - _binary_default_settings_yml_start);
}

// Bail out on errors.
static void error(const char *message)
{
    std::cerr << "CONFIGURATION ERROR: " << message << "\n";
    exit(1);
}

// Conditionally assign a variable if the node is defined and scalar.
template <typename T>
static bool assign(T &var, const YAML::Node &node)
{
    if (node.IsScalar()) {
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

// Stringify a boolean
static const char * str(bool x)
{
    return x ? "true" : "false";
}

// Stringify a vector
static const std::string str(Pt::pt3D v)
{
    return std::string("[")
        + std::to_string(v.x()) + ", "
        + std::to_string(v.y()) + ", "
        + std::to_string(v.z()) + "]";
}


/***********************************************************************
 * Public API.
 */

Settings::Settings()
{
    verbose = false;
    withTsv = true;
    read(YAML::Load(get_default_yaml()));  // load defaults
}

void Settings::dumpDefaults()
{
    std::cout << get_default_yaml();
}

void Settings::infos()
{
    std::cout << "outputs:\n";
    std::cout << "  directory: " << r_path_output_dir << "\n";
    std::cout << "  file_basename: " << simName << "\n";
    std::cout << "  evol_time_step: " << time_step << "\n";
    std::cout << "  evol_columns:\n";
    for (auto it = evol_columns.begin(); it != evol_columns.end(); ++it) {
        std::cout << "    - " << *it << "\n";
    }
    std::cout << "  evol_header: " << str(evol_header) << "\n";
    std::cout << "  take_photo: " << save_period << "\n";
    std::cout << "  vtk_file: " << str(withVtk) << "\n";
    std::cout << "mesh:\n";
    std::cout << "  filename: " << pbName << "\n";
    std::cout << "  scaling_factor: " << _scale << "\n";
    std::cout << "  volume_regions:\n";
    for (auto it = paramTetra.begin(); it != paramTetra.end(); ++it) {
        if (it->regName == "__default__")  // skip
            continue;
        std::cout << "    " << it->regName << ":\n";
        std::cout << "      Ae: " << it->A << "\n";
        std::cout << "      Js: " << it->J << "\n";
        std::cout << "      K: " << it->K << "\n";
        if (it->K != 0)
            std::cout << "      uk: " << str(it->uk) << "\n";
        std::cout << "      K3: " << it->K3 << "\n";
        if (it->K3 != 0) {
            std::cout << "      ex: " << str(it->ex) << "\n";
            std::cout << "      ey: " << str(it->ey) << "\n";
            std::cout << "      ez: " << str(it->ez) << "\n";
        }
        std::cout << "      alpha_LLG: " << it->alpha_LLG << "\n";
    }
    std::cout << "  surface_regions:\n";
    for (auto it = paramFacette.begin(); it != paramFacette.end(); ++it) {
        if (it->regName == "__default__")  // skip
            continue;
        std::cout << "    " << it->regName << ":\n";
        std::cout << "      suppress_charges: "
            << str(it->suppress_charges) << "\n";
        std::cout << "      Ks: " << it->Ks << "\n";
        if (it->Ks != 0)
            std::cout << "      uk: " << str(it->uk) << "\n";
    }
    std::cout << "initial_magnetization: ";
    if (restoreFileName.empty())
        std::cout << "[\""
            << sMx << "\", \"" << sMy << "\", \"" << sMz << "\"]\n";
    else
        std::cout << restoreFileName << "\n";
    std::cout << "recentering:\n";
    std::cout << "  enable: " << str(recenter) << "\n";
    if (recenter) {
        std::cout << "  direction: ";
        switch (recentering_direction) {
            case Pt::IDX_UNDEF: std::cout << "UNDEF\n"; break; 
            case Pt::IDX_X: std::cout << "X\n"; break;
            case Pt::IDX_Y: std::cout << "Y\n"; break;
            case Pt::IDX_Z: std::cout << "Z\n"; break;
        }
        std::cout << "  threshold: " << threshold << "\n";
    }
    std::cout << "Bext: [\""
        << sBx << "\", \"" << sBy << "\", \"" << sBz << "\"]\n";
    std::cout << "spin_transfer_torque:\n";
    std::cout << "  enable: " << str(stt_flag) << "\n";
    if (stt_flag) {
        std::cout << "  sigma: " << p_stt.sigma << "\n";
        std::cout << "  dens_state: " << p_stt.N0 << "\n";
        std::cout << "  beta: " << p_stt.beta << "\n";
        std::cout << "  l_J: " << p_stt.lJ << "\n";
        std::cout << "  l_sf: " << p_stt.lsf << "\n";
        std::cout << "  V_file: " << str(p_stt.V_file) << "\n";
        std::cout << "  boundary_conditions:";
        if (p_stt.boundaryCond.size() == 0) {
            std::cout << " {}\n";  // empty map
        } else {
            std::cout << "\n";
            for (unsigned int i = 0; i < p_stt.boundaryCond.size(); i++) {
                std::cout << "    \"" << p_stt.boundaryCond[i].first
                    << "\": " << p_stt.boundaryCond[i].second << "\n";
            }
        }
    }
    
    std::cout << "finite_element_solver:\n";
    std::cout << "  nb_threads: " << solverNbTh << "\n";
    std::cout << "  max(iter): " << MAXITER << "\n";
    std::cout << "  refresh_preconditioner_every: " << REFRESH_PRC << "\n";
    std::cout << "demagnetizing_field_solver:\n";
    std::cout << "  nb_threads: " << scalfmmNbTh << "\n";
    std::cout << "time_integration:\n";
    std::cout << "  final_time: " << tf << "\n";
    std::cout << "  max(du): " << DUMAX << "\n";
    std::cout << "  min(dt): " << dt_min << "\n";
    std::cout << "  max(dt): " << dt_max << "\n";
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
        }
        assign(simName, outputs["file_basename"]);
        assign(withVtk, outputs["vtk_file"]);
        assign(time_step, outputs["evol_time_step"]);
        YAML::Node take_photo = outputs["take_photo"];
        if (take_photo.Scalar() == "true") {  // catch an easy mistake
            error("outputs.take_photo should be an integer or `false'.");
        } else if (take_photo.Scalar() == "false") {
            save_period = 0;
        } else {
            if (assign(save_period, take_photo) && save_period < 0)
                { save_period = 0; }
        }
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
            int default_idx = findTetraRegionIdx("__default__");
            Tetra::prm *default_volume = nullptr;
            if (default_idx > -1)
                default_volume = &paramTetra[default_idx];
            for (auto it = volumes.begin(); it != volumes.end(); ++it) {
                std::string name = it->first.as<std::string>();
                YAML::Node volume = it->second;
                Tetra::prm p;
                if (default_volume) p = *default_volume;
                p.regName = name;
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
                
                paramTetra.push_back(p);
            }
        }  // mesh.volume_regions
        YAML::Node surfaces = mesh["surface_regions"];
        if (surfaces) {
            if (!surfaces.IsMap())
                error("mesh.surface_regions should be a map.");
            int default_idx = findFacetteRegionIdx("__default__");
            Facette::prm *default_surface = nullptr;
            if (default_idx > -1)
                default_surface = &paramFacette[default_idx];
            for (auto it = surfaces.begin(); it != surfaces.end(); ++it) {
                std::string name = it->first.as<std::string>();
                YAML::Node surface = it->second;
                Facette::prm p;
                if (default_surface) p = *default_surface;
                p.regName = name;
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
        assign(recenter, recentering["enable"]);
        std::string direction = recentering["direction"].as<std::string>("Z");
        if (direction.length() != 1 || direction[0] < 'X' || direction[0] > 'Z')
            error("recentering.direction should be X, Y or Z.");
        switch (direction[0]) {
            case 'X': recentering_direction = Pt::IDX_X; break;
            case 'Y': recentering_direction = Pt::IDX_Y; break;
            case 'Z': recentering_direction = Pt::IDX_Z; break;
        }
        assign(threshold, recentering["threshold"]);
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
        assign(stt_flag, stt["enable"]);
        assign(p_stt.sigma, stt["sigma"]);
        assign(p_stt.N0, stt["dens_state"]);
        assign(p_stt.beta, stt["beta"]);
        assign(p_stt.lJ, stt["l_J"]);
        assign(p_stt.lsf, stt["l_sf"]);
        assign(p_stt.V_file, stt["V_file"]);
        
        YAML::Node bound_cond = stt["boundary_conditions"];
        if (bound_cond)
        	{
        	if (!bound_cond.IsMap())
                	error("stt.boundary_conditions should be a map.");
        	else
        		{
        		for (auto it = bound_cond.begin(); it != bound_cond.end(); ++it)
        			{
                		std::string name = it->first.as<std::string>();
                		double V = it->second.as<double>();
        			p_stt.boundaryCond.push_back(std::make_pair(name,V)); 
        			//std::cout<< "surf name : " << name << "\tV = " << V <<std::endl;
        			}
        		}
        	}
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
        assign(tf, time_integration["final_time"]);
    }  // time_integration

    // outputs.file_basename defaults to base name of mesh.filename.
    if (simName.empty() && !pbName.empty()) {
        simName = pbName;

        // Remove directory part.
        size_t pos = simName.rfind('/');
        if (pos != simName.npos)
            simName.erase(0, pos + 1);

        // Remove everything after last dot.
        pos = simName.rfind('.');
        if (pos != simName.npos)
            simName.erase(pos);
    }
}

bool Settings::read(std::string filename)
{
    YAML::Node config;
    if (filename == "-")
        config = YAML::Load(std::cin);
    else
        config = YAML::LoadFile(filename);
    if (config.IsNull())
        return false;
    read(config);
    return true;
}
