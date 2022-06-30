#include <iostream>
#include <fstream>
#include <string>
#include <ctime>
#include <cctype>
#include <exception>
#include <stdio.h>  // for perror()
#include <unistd.h>  // for sysconf(), stat()
#include <sys/stat.h>  // for mkdir(), stat()
#include <sys/types.h>  // for mkdir(), stat()

/* Silence a warning internal to Boost, fixed in Boost 1.76.0. */
#if BOOST_VERSION < 107600
# define BOOST_BIND_GLOBAL_PLACEHOLDERS
#endif

#include<boost/property_tree/ptree.hpp>
#include <boost/property_tree/json_parser.hpp>

#include "feellgoodSettings.h"

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


void Settings::read(std::string fileJson)
{
boost::property_tree::ptree root;
boost::property_tree::ptree sub_tree,s_sub_tree;
    
boost::property_tree::read_json(fileJson,root);
	
try { sub_tree = root.get_child("outputs"); }
    catch (std::exception &e)
        { std::cout << e.what() << std::endl; }

r_path_output_dir = sub_tree.get<std::string>("directory", ".");
// Normalize directory name.
if (r_path_output_dir.empty())
    { r_path_output_dir = "."; }
if (r_path_output_dir.length() > 1 && r_path_output_dir.back() == '/')
    { r_path_output_dir.pop_back(); }

create_dir_if_needed(r_path_output_dir);

simName = sub_tree.get<std::string>("file_basename");
withVtk = sub_tree.get<bool>("vtk_file",0);
time_step = sub_tree.get<double>("evol_time_step",1e-7);
save_period = sub_tree.get<int>("take_photo",0);

try {s_sub_tree = sub_tree.get_child("evol_columns");}
catch(std::exception &e)
    { std::cout << e.what() << std::endl; }

evol_columns.clear();
for(boost::property_tree::ptree::value_type &cell :s_sub_tree)
    { evol_columns.push_back( cell.second.get_value<std::string>() ); }
    
evol_header = sub_tree.get<bool>("evol_header",false);

try { sub_tree = root.get_child("mesh"); }
    catch (std::exception &e)
        { std::cout << e.what() << std::endl; }

pbName = sub_tree.get<std::string>("filename");
double s = sub_tree.get<double>("scaling_factor",0.0);

if (s <= 0.0)
    {
    std::cerr << "scaling factor must be defined and strictly positive in json settings file" << std::endl;
    exit(1);
    }
setScale(s);
    
try { s_sub_tree = sub_tree.get_child("volume_regions"); }
catch (std::exception &e)
    { std::cout << e.what() << std::endl; }
	
for (boost::property_tree::ptree::value_type &s : s_sub_tree)
    {
    std::string name_reg = s.first;
    if (!name_reg.empty())
        {
        Tetra::prm p;
        p.reg = stoi(name_reg);
        for (boost::property_tree::ptree::value_type &sub_k : s_sub_tree.get_child(name_reg))
            {
            if (sub_k.first == "Ae") {p.A = sub_k.second.get_value<double>();}
            if (sub_k.first == "alpha_LLG") {p.alpha_LLG = sub_k.second.get_value<double>();}
            
            if (sub_k.first == "K")
                { p.K = sub_k.second.get_value<double>(); }
            
            if (sub_k.first == "uk")
                { p.uk = readUnitVector(sub_k,"uk"); }
            
            if (sub_k.first == "K3")
                { p.K3 = sub_k.second.get_value<double>(); }
            
            if (sub_k.first == "ex")
                { p.ex = readUnitVector(sub_k,"ex"); }
                
            if (sub_k.first == "ey")
                { p.ey = readUnitVector(sub_k,"ey"); }
            
            if (sub_k.first == "ez")
                { p.ez = readUnitVector(sub_k,"ez"); }
            
            if ( !Pt::isOrthogonal(p.ex,p.ey,p.ez,USER_TOL) )
                { std::cout << "Warning : {ex,ey,ez} is not orthogonal" << std::endl; }
            
            if (sub_k.first == "Js")
                { p.J = sub_k.second.get_value<double>(); }
            }
        p.p_STT.beta = 0.0;
        p.p_STT.N0 = 1.0;
        p.p_STT.sigma = 1.0;
        p.p_STT.lJ = 1.0;
        p.p_STT.lsf = 1.0;
        p.p_STT.gamma0 = 1.0;
        p.p_STT.func = [](Pt::pt3D) {return 1;};
        
        paramTetra.push_back(p);
        }
    }
		
try { s_sub_tree = sub_tree.get_child("surface_regions"); }
catch (std::exception &e)
    { std::cout << e.what() << std::endl; }

for (boost::property_tree::ptree::value_type &s : s_sub_tree)
    {
    std::string name_reg = s.first;
    if (!name_reg.empty())
        {
        Facette::prm p;
        p.reg = stoi(name_reg);
        for (boost::property_tree::ptree::value_type &sub_k : s_sub_tree.get_child(name_reg))
            {
            if (sub_k.first == "Ks")
                { p.Ks = sub_k.second.get_value<double>(); }
            if (sub_k.first == "uk")
                { p.uk = readUnitVector(sub_k,"uk"); }
            if (sub_k.first == "suppress_charges") 
                { p.suppress_charges = sub_k.second.get_value<bool>(); }
            }
        paramFacette.push_back(p);
        //p.infos();	
        }
    }

restoreFileName = root.get("initial_magnetization","");
if(restoreFileName == "")
    {
        try { sub_tree = root.get_child("initial_magnetization"); }
            catch (std::exception &e)
                { std::cout << e.what() << std::endl; }

        if (sub_tree.size() != 3)
            {
            std::cerr << "initial_magnetization should have three components.\n";
            exit(1);
            }
        auto it = sub_tree.begin();
        sMx = it->second.data(); it++;
        sMy = it->second.data(); it++;
        sMz = it->second.data();
        doCompile3Dprm();
        
        if (verbose) { std::cout << "initial magnetization defined from a math expression" << std::endl; }
    }
else if (verbose) { std::cout<< "initial magnetization defined from file :" << restoreFileName <<std::endl; }    

try 
    {
    sub_tree = root.get_child("recentering"); 
    recenter = true;
    char recentering_dir = sub_tree.get<char>("direction",'Z');

    if((recentering_dir != 'X') && (recentering_dir != 'Y') && (recentering_dir != 'Z'))
        { std::cout << "unknown recentering direction !"<< std::endl; exit(1); }
    if (recentering_dir == 'X') { recentering_direction = Pt::IDX_X; }
        else if (recentering_dir == 'Y') { recentering_direction = Pt::IDX_Y; }
            else { recentering_direction = Pt::IDX_Z; }

    threshold = sub_tree.get<double>("threshold",0.1);
    }
    catch (std::exception &e)
        { //std::cout << e.what() << std::endl; 
        recenter = false;
        }

try {sub_tree = root.get_child("Bext");}
catch(std::exception &e)
    { std::cout << e.what() << std::endl; }

if (sub_tree.size() != 3)
    {
    std::cerr << "Bext should have three components.\n";
    exit(1);
    }
auto it = sub_tree.begin();
sBx = it->second.data(); it++;
sBy = it->second.data(); it++;
sBz = it->second.data();
doCompile1Dprm();

try
    {
    sub_tree = root.get_child("spin_transfer_torque");
    
    p_stt.func = [](Pt::pt3D) {return 1;};
    
    try
        {
        p_stt.gamma0 = sub_tree.get<double>("gamma0",0.0);
        p_stt.sigma = sub_tree.get<double>("sigma",1.0);
        p_stt.N0 = sub_tree.get<double>("dens_state",1.0);
        p_stt.beta = sub_tree.get<double>("beta",0.0);
    
        p_stt.lJ = sub_tree.get<double>("l_J",1.0);
        p_stt.lsf = sub_tree.get<double>("l_sf",1.0);
        std::string reg_str = sub_tree.get<std::string>("volume_region_reference"); 
        p_stt.reg = stoi(reg_str);
        p_stt.bc = Tetra::boundary_conditions::Undef; // default value
        try 
            {
            s_sub_tree = sub_tree.get_child("boundary_conditions");
            
            for (boost::property_tree::ptree::value_type &s : s_sub_tree)
                {
                std::string name_reg = s.first;
                if (!name_reg.empty())
                    {
                    Tetra::bc_data bc_vals;
                    bc_vals.reg = stoi(name_reg);
                    
                    try { bc_vals.val = s_sub_tree.get_child(s.first).get<double>("V"); bc_vals.typ = Tetra::type_val_reg::potV; }
                    catch(std::exception &e) 
                        {  
                        try { bc_vals.val = s_sub_tree.get_child(s.first).get<double>("J"); bc_vals.typ = Tetra::type_val_reg::densJ; }
                        catch(std::exception &e) { std::cout << "unrecognized boundary condition physical constant on region " << bc_vals.reg << std::endl; exit(1); }
                        }
                    
                    p_stt.full_bc_data.push_back( bc_vals );
                    }
                }   
            }
        catch(std::exception &e)
            { std::cout << "warning : missing boundary conditions in spin_transfer_torque" << std::endl; }
        p_stt.define_boundary_conditions();
        if (p_stt.bc != Tetra::boundary_conditions::Undef ) {stt_flag = true;}
        }
    catch(std::exception &e)
        { stt_flag = false; std::cout<< "incorrect syntax/missing parameter in spin_transfer_torque section" <<std::endl; }
    }
catch(std::exception &e)
    { stt_flag = false; }

// The number of available processors (actually, hardware threads) is
// the default for the number of threads to spin.
int available_cpu_count = sysconf(_SC_NPROCESSORS_ONLN);

try { sub_tree = root.get_child("demagnetizing_field_solver"); }
catch (std::exception &e)
    { std::cout << e.what() << std::endl; }
scalfmmNbTh = sub_tree.get<int>("nb_threads",0);
if (scalfmmNbTh <= 0)
    { scalfmmNbTh = available_cpu_count; }

try { sub_tree = root.get_child("finite_element_solver"); }
catch (std::exception &e)
    { std::cout << e.what() << std::endl; }
solverNbTh = sub_tree.get<int>("nb_threads",0);
if (solverNbTh <= 0)
    { solverNbTh = available_cpu_count; }

MAXITER = sub_tree.get<int>("max(iter)",500);
REFRESH_PRC = sub_tree.get<int>("refresh_preconditioner_every",20);

try { sub_tree = root.get_child("time_integration"); }
catch (std::exception &e)
    { std::cout << e.what() << std::endl; }
    
DUMAX = sub_tree.get<double>("max(du)",0.02);
dt_min = sub_tree.get<double>("min(dt)",1e-14);
dt_max = sub_tree.get<double>("max(dt)",1e-7);
tf = sub_tree.get<double>("final_time",dt_max);
}

Pt::pt3D Settings::readUnitVector(boost::property_tree::ptree::value_type &sub_k, const std::string varName)
{
Pt::pt3D vect;

std::vector<double> X;
for(boost::property_tree::ptree::value_type &row : sub_k.second) { X.push_back( row.second.get_value<double>() ); }
if(X.size() != Pt::DIM) {std::cout<<"wrong number of components for " << varName <<std::endl; exit(1);}
vect = Pt::pt3D(X[0],X[1],X[2]);
if ( fabs(vect.norm()-1.0) > USER_TOL )
    {std::cout<<"WARNING: " << varName <<"= " << vect << " not a unit vector"<<std::endl; }
vect.normalize();
return vect;
}
