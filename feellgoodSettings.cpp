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
std::cout << "\t final time of the simulation\t\t" << tf << std::endl;
std::cout << "\t initial time step\t" << dt << std::endl << std::endl;
std::cout << "\t save energy values every " << n1 << " iterations" << std::endl;
std::cout << "\t snapshot of the magnetization configuration every " << n2 << " iterations" << std::endl;
if (restore)
    { std::cout << "\t restore initial magnetization distribution from file : " << restoreFileName <<std::endl; }
else
    { std::cout << "\t initial magnetization distribution from math expression" << std::endl; }

std::cout << "\t applied field Bext = [ " << Bext[0] << ",\t" << Bext[1] << ",\t" << Bext[2] << " ] T" << std::endl;

for(unsigned int i=0;i<paramTetra.size();i++) {paramTetra[i].infos();}
for(unsigned int i=0;i<paramFacette.size();i++) {paramFacette[i].infos();}

std::cout << solverNbTh+1 << " threads for assembling matrix." << std::endl;
}


void Settings::read(std::string fileJson)
{
boost::property_tree::ptree root;
boost::property_tree::ptree sub_tree,s_sub_tree;
    
boost::property_tree::read_json(fileJson,root);
	
std::cout << "\nparsing "<< fileJson <<" :" << std::endl;
	
    
try { sub_tree = root.get_child("outputs"); }
    catch (std::exception &e)
        { std::cout << e.what() << std::endl; }

r_path_output_dir = sub_tree.get<std::string>("directory");
simName = sub_tree.get<std::string>("file basename");
withVtk = sub_tree.get<bool>("vtk file",0);
n1 = sub_tree.get<int>("save_energies",0);
n2 = sub_tree.get<int>("take_photo",0);
verbose = sub_tree.get<bool>("verbose",true);

try {s_sub_tree = sub_tree.get_child("evol columns");}
catch(std::exception &e)
    { std::cout << e.what() << std::endl; }

evol_columns.clear();
for(boost::property_tree::ptree::value_type &cell :s_sub_tree)
    { evol_columns.push_back( cell.second.get_value<std::string>() ); }
    
try { sub_tree = root.get_child("mesh"); }
    catch (std::exception &e)
        { std::cout << e.what() << std::endl; }

pbName = sub_tree.get<std::string>("filename");
double s = sub_tree.get<double>("scaling factor",0.0);
EPSILON = sub_tree.get<double>("epsilon",1e-16);

if (s <= 0.0)
    {
    std::cerr << "scaling factor must be defined and strictly positive in json settings file" << std::endl;
    exit(1);
    }
setScale(s);
    
std::cout << "\tvolumic regions..." << std::endl;


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
            if (sub_k.first == "alpha") {p.alpha = sub_k.second.get_value<double>();}
            
            if (sub_k.first == "Ka")
                {
                p.K = sub_k.second.get_value<double>();
                if (p.K == 0) { std::cout<< "\tKa is zero, flag NO_ANISOTROPY to true" << std::endl;}
                }
            if (sub_k.first == "a")
                {int i=0;
                for(boost::property_tree::ptree::value_type &row : sub_k.second)
                    {
                    int j=0;
                    for(boost::property_tree::ptree::value_type &coeff : row.second)
                        {
                        p.uk[i][j] = coeff.second.get_value<double>();	
                        j++;
                        }
                    i++;	
                    }
                }
			if (sub_k.first == "Js") {p.J = sub_k.second.get_value<double>();}
            }
        paramTetra.push_back(p);
        //p.infos();	
        }
    }
		
std::cout << "\tsurface regions..." << std::endl;

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
                if (p.Ks == 0) { std::cout<< "\tKs is zero, flag NO_SURF_ANISOTROPY to true" << std::endl;}
                }
            if (sub_k.first == "uk")
                {int i=0;
                for(boost::property_tree::ptree::value_type &row : sub_k.second)
                    {
                    p.uk[i] = row.second.get_value<double>();
                    i++;	
                    }
                }
			if (sub_k.first == "Js") {p.Js = sub_k.second.get_value<double>();}	
            }
        paramFacette.push_back(p);
        //p.infos();	
        }
    }

restore = root.get<bool>("restore",0);
if (restore)
    { restoreFileName = root.get<std::string>("restore from file"); }
else
    {
    try { sub_tree = root.get_child("initial magnetization"); }
    catch (std::exception &e)
        { std::cout << e.what() << std::endl; }

    sMx = sub_tree.get<std::string>("Mx");
    sMy = sub_tree.get<std::string>("My");
    sMz = sub_tree.get<std::string>("Mz");
    doCompile();
    }

try { sub_tree = root.get_child("recentering"); }
catch (std::exception &e)
        { std::cout << e.what() << std::endl; }

recenter = sub_tree.get<bool>("recenter",false);
recentering_direction = sub_tree.get<char>("direction",'Z');
threshold = sub_tree.get<double>("threshold",0.1);

try {sub_tree = root.get_child("Bext");}
catch(std::exception &e)
    { std::cout << e.what() << std::endl; }
int j=0;    
for(boost::property_tree::ptree::value_type &cell :sub_tree)
    {
    Bext[j] = cell.second.get_value<double>();
    j++;
    }

try {sub_tree = root.get_child("spin polarized current");}
catch(std::exception &e)
    { std::cout << e.what() << std::endl; }
Uz = sub_tree.get<double>("Uz",0.0);
beta = sub_tree.get<double>("beta",0.0);

try { sub_tree = root.get_child("demagnetization field solver"); }
catch (std::exception &e)
    { std::cout << e.what() << std::endl; }
scalfmmNbTh = sub_tree.get<int>("nbThreads",8);
analytic_corr = sub_tree.get<bool>("analytic corrections",true);

try { sub_tree = root.get_child("finite element solver"); }
catch (std::exception &e)
    { std::cout << e.what() << std::endl; }
solverNbTh = sub_tree.get<int>("nbThreads",8);

MAXITER = sub_tree.get<int>("max(iter)",500);
REFRESH_PRC = sub_tree.get<int>("refresh preconditionner every",20);

try { sub_tree = root.get_child("time integration"); }
catch (std::exception &e)
    { std::cout << e.what() << std::endl; }
theta = sub_tree.get<double>("theta",0.5);
DUMIN = sub_tree.get<double>("min(du)",1e-9);
DUMAX = sub_tree.get<double>("max(du)",0.02);
DTMIN = sub_tree.get<double>("min(dt)",1e-14);
DTMAX = sub_tree.get<double>("max(dt)",1e-7);
TAUR = 100.*DTMAX; // pourquoi 100 ?
dt = sub_tree.get<double>("initial dt",1e-9);
tf = sub_tree.get<double>("final_time",0);
}
