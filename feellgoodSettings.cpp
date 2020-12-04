#include <iostream>
#include <fstream>
#include <string>
#include <ctime>
#include <cctype>
#include <exception>

#include<boost/property_tree/ptree.hpp>
#include <boost/property_tree/json_parser.hpp>

#include "feellgoodSettings.h"

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

std::cout << "recentering along " << recentering_direction << " with threshold = " << threshold << std::endl;

for(unsigned int i=0;i<paramTetra.size();i++) {paramTetra[i].infos();}
for(unsigned int i=0;i<paramFacette.size();i++) {paramFacette[i].infos();}

std::cout << solverNbTh+1 << " threads for assembling matrix." << std::endl;
}


timing Settings::read(std::string fileJson)
{
boost::property_tree::ptree root;
boost::property_tree::ptree sub_tree,s_sub_tree;
    
boost::property_tree::read_json(fileJson,root);
	
try { sub_tree = root.get_child("outputs"); }
    catch (std::exception &e)
        { std::cout << e.what() << std::endl; }

r_path_output_dir = sub_tree.get<std::string>("directory");
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
                {
                double X[Pt::DIM];
                int i=0;
                for(boost::property_tree::ptree::value_type &row : sub_k.second) { X[i] = row.second.get_value<double>();i++; }
                p.uk = Pt::pt3D(X[0],X[1],X[2]);
                }
            
            if (sub_k.first == "K3")
                { p.K3 = sub_k.second.get_value<double>(); }
            
            if (sub_k.first == "alpha")
                {
                double X[Pt::DIM];
                int i=0;
                for(boost::property_tree::ptree::value_type &row : sub_k.second) { X[i] = row.second.get_value<double>();i++; }
                p.alpha = Pt::pt3D(X[0],X[1],X[2]);
                }
            
            if (sub_k.first == "beta")
                {
                double X[Pt::DIM];
                int i=0;
                for(boost::property_tree::ptree::value_type &row : sub_k.second) { X[i] = row.second.get_value<double>();i++; }
                p.beta = Pt::pt3D(X[0],X[1],X[2]);
                }
            
            if (sub_k.first == "gamma")
                {
                double X[Pt::DIM];
                int i=0;
                for(boost::property_tree::ptree::value_type &row : sub_k.second) { X[i] = row.second.get_value<double>();i++; }
                p.gamma = Pt::pt3D(X[0],X[1],X[2]);
                }
            
            if (sub_k.first == "Js") {p.J = sub_k.second.get_value<double>();}
            }
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
                {
                p.Ks = sub_k.second.get_value<double>();
                }
            if (sub_k.first == "uk")
                {double X[Pt::DIM];
		int i=0;
                for(boost::property_tree::ptree::value_type &row : sub_k.second)
                    { X[i] = row.second.get_value<double>();i++; }
		p.uk = Pt::pt3D(X[0],X[1],X[2]);                
		}
	if (sub_k.first == "suppress_charges") {p.suppress_charges = sub_k.second.get_value<bool>();}
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

        sMx = sub_tree.get<std::string>("Mx");
        sMy = sub_tree.get<std::string>("My");
        sMz = sub_tree.get<std::string>("Mz");
        doCompile3Dprm();
        
        if (verbose) { std::cout << "initial magnetization defined from a math expression" << std::endl; }
    }
else if (verbose) { std::cout<< "initial magnetization defined from file :" << restoreFileName <<std::endl; }    

try 
    {
    sub_tree = root.get_child("recentering"); 
    recentering_direction = sub_tree.get<char>("direction",'Z');

    if((recentering_direction != 'X') && (recentering_direction != 'Y') && (recentering_direction != 'Z'))
        { std::cout << "unknown recentering direction !"<< std::endl; exit(1); }
    threshold = sub_tree.get<double>("threshold",0.1);
    }
    catch (std::exception &e)
        { //std::cout << e.what() << std::endl; 
        recenter = false;
        }

try {sub_tree = root.get_child("Bext");}
catch(std::exception &e)
    { std::cout << e.what() << std::endl; }

sBx = sub_tree.get<std::string>("Bx");
sBy = sub_tree.get<std::string>("By");
sBz = sub_tree.get<std::string>("Bz");
doCompile1Dprm();
    
/*
std::vector<double> val_vect;
for(boost::property_tree::ptree::value_type &cell :sub_tree)
    { val_vect.push_back( nu0 * cell.second.get_value<double>() ); }

if(val_vect.size() != Pt::DIM) {std::cout<<"wrong number of field components"<<std::endl;}
else
    { Hext = Pt::pt3D(val_vect[0],val_vect[1],val_vect[2]); }
*/
    
try
    {
    sub_tree = root.get_child("spin_polarized_current");
    Uz = sub_tree.get<double>("Uz",0.0);
    beta = sub_tree.get<double>("beta",0.0);
    }
catch(std::exception &e)
    { std::cout << "no spin polarized" << std::endl; }


try { sub_tree = root.get_child("demagnetization_field_solver"); }
catch (std::exception &e)
    { std::cout << e.what() << std::endl; }
scalfmmNbTh = sub_tree.get<int>("nb_threads",8);

try { sub_tree = root.get_child("finite_element_solver"); }
catch (std::exception &e)
    { std::cout << e.what() << std::endl; }
solverNbTh = sub_tree.get<int>("nb_threads",8);

MAXITER = sub_tree.get<int>("max(iter)",500);
REFRESH_PRC = sub_tree.get<int>("refresh_preconditionner_every",20);

try { sub_tree = root.get_child("time_integration"); }
catch (std::exception &e)
    { std::cout << e.what() << std::endl; }
    
DUMIN = sub_tree.get<double>("min(du)",1e-9);
DUMAX = sub_tree.get<double>("max(du)",0.02);

double dt_min = sub_tree.get<double>("min(dt)",1e-14);
double dt_max = sub_tree.get<double>("max(dt)",1e-7);
double tf = sub_tree.get<double>("final_time",dt_max);

timing t_prm = timing(0.0,tf,dt_min,dt_max);

return(t_prm);
}
