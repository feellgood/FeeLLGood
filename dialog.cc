#include <iostream>
#include <fstream>
#include <string>
#include <ctime>
#include <cctype>
#include <unistd.h> // for getpid()

#include <exception>

#include<boost/property_tree/ptree.hpp>
#include <boost/property_tree/json_parser.hpp>

#include "feellgoodSettings.h"

void Settings::helloDialog()
{
std::cout << "\n\t ******************************\n";
std::cout <<   "\t *         feeLLGood          *\n";
std::cout <<   "\t *        version 2018        *\n";
std::cout <<   "\t *      cnrs Grenoble-INP     *\n";
std::cout <<   "\t ******************************\n";
std::cout << "\t process\t\t" << getpid() << std::endl;
}

void Settings::printToTerminal(std::vector<Seq> &seq)
{
std::cout << "\t name simul : "<< simName << std::endl;
std::cout << "\t name mesh file : " << pbName << std::endl;
std::cout << "\t temps final\t\t" << tf << std::endl;
std::cout << "\t pas de temps init\t" << dt << std::endl << std::endl;
std::cout << "\t sauve energies chaque " << n1 << " iterations" << std::endl;
std::cout << "\t photo config chaque " << n2 << " iterations" << std::endl;
std::cout << "\t restauration a partir d'un fichier sol.in: " << restore << std::endl <<std::endl;
std::cout << "\t sequences de champ applique Bext" << std::endl;
for (std::vector<Seq>::iterator it = seq.begin(); it!=seq.end(); ++it) {
    std::cout << it->Bini << "\t" << it->Bfin << "\t" << it->dB << "\t" << it->a[0] << "\t" << it->a[1] << "\t" << it->a[2] << std::endl;
    }
for(unsigned int i=0;i<paramTetra.size();i++) {paramTetra[i].infos();}
for(unsigned int i=0;i<paramFacette.size();i++) {paramFacette[i].infos();}
}

void Settings::extract_comment(std::istream &flux)
{
std::string comment;
char c;
while (1)
      {
      while (flux.get(c))  // lecture des caracteres blancs (espaces, nl, tab, cr, np)
	    if (isspace(c)==0)
	       {
	       flux.putback(c);
	       break;
	       }

      flux.get(c);  
      if (c=='#') 
	 getline(flux, comment,'\n');  // lecture du commentaire
      else
	 {
	 flux.putback(c);
	 break;
	 }
      }
}

void Settings::dialog(std::vector<Seq> &seq)
{
int SEQ;

std::cout << "\t name simul ? "<< std::endl;
extract_comment(std::cin);   std::cin >> simName;
std::cout << "\t name simul\t\t"<< simName << std::endl;

std::cout << "\t file mesh ?"<< std::endl;
extract_comment(std::cin);   std::cin >> pbName;
std::cout << "\t name file pro\t" << pbName << std::endl;

std::ifstream fin(pbName);
if (!fin){
    std::cerr << "Cannot open problem file, fileName = " << pbName << std::endl;
    exit(1);}
fin.close();

extract_comment(std::cin);              std::cin >> tf;
extract_comment(std::cin);              std::cin >> dt;
extract_comment(std::cin);              std::cin >> theta;
extract_comment(std::cin);              std::cin >> n1;
extract_comment(std::cin);              std::cin >> n2;
extract_comment(std::cin);              std::cin >> restore;
extract_comment(std::cin);              std::cin >> SEQ;

for (int nseq=0; nseq<SEQ; nseq++)
	{
    Seq s;
    extract_comment(std::cin);          std::cin >> s.Bini >> s.Bfin >> s.dB >> s.a[0] >> s.a[1] >> s.a[2];
    s.dB=fabs(s.dB)*(s.Bfin<s.Bini? -1.0: 1.0);
    normalize(s.a);
    seq.push_back(s);
    }
printToTerminal(seq);
}

void Settings::read(std::string fileJson,std::vector<Seq> &seq)
{
	boost::property_tree::ptree root;
	boost::property_tree::read_json(fileJson,root);
	
	std::cout << "parsing parameters and settings from json file :\n" << std::endl;
	
	r_path_output_dir = root.get<std::string>("output directory");
	simName = root.get<std::string>("output file basename");
	pbName = root.get<std::string>("mesh filename");
	tf = root.get<double>("final_time",0);
	
	n1 = root.get<int>("save_energies",0);
	n2 = root.get<int>("take_photo",0);
	restore = root.get<bool>("restore",0);
	restoreFileName = root.get<std::string>("restore from file");
    recentering = root.get<bool>("recentering",true);
    double trucs[6];
	for (boost::property_tree::ptree::value_type &s : root.get_child("field_sequence"))
		{
		int j=0;
		for(boost::property_tree::ptree::value_type &cell : s.second)
			{
			trucs[j] = cell.second.get_value<double>();
			j++;
			}
		Seq field;
		field.Bini = trucs[0];
		field.Bfin = trucs[1];
		field.dB = trucs[2];
		field.a[0] = trucs[3];
		field.a[1] = trucs[4];
		field.a[2] = trucs[5];
		normalize(field.a);
		seq.push_back(field);
		}
std::cout << "volumic regions..." << std::endl;

boost::property_tree::ptree sub_tree;

try { sub_tree = root.get_child("volume_regions"); }
catch (std::exception &e)
    { std::cout << e.what() << std::endl; }
	
	for (boost::property_tree::ptree::value_type &s : sub_tree)
		{
		std::string name_reg = s.first;
		if (!name_reg.empty())
			{
			Tetra::prm p;
			p.reg = stoi(name_reg);
			for (boost::property_tree::ptree::value_type &sub_k : sub_tree.get_child(name_reg))
				{
				if (sub_k.first == "Ae") {p.A = sub_k.second.get_value<double>();}
				if (sub_k.first == "alpha") {p.alpha = sub_k.second.get_value<double>();}
				
				if (sub_k.first == "Ka")
					{
					p.K = sub_k.second.get_value<double>();
					if (p.K == 0) { std::cout<< "Ka is zero, flag NO_ANISOTROPY to true" << std::endl;}
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
		
std::cout << "surfacic regions..." << std::endl;

try { sub_tree = root.get_child("surface_regions"); }
catch (std::exception &e)
    { std::cout << e.what() << std::endl; }

for (boost::property_tree::ptree::value_type &s : sub_tree)
		{
		std::string name_reg = s.first;
		if (!name_reg.empty())
			{
			Facette::prm p;
			p.reg = stoi(name_reg);
			for (boost::property_tree::ptree::value_type &sub_k : sub_tree.get_child(name_reg))
				{
				
				if (sub_k.first == "Ks")
					{
					p.Ks = sub_k.second.get_value<double>();
					if (p.Ks == 0) { std::cout<< "Ks is zero, flag NO_SURF_ANISOTROPY to true" << std::endl;}
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

try { sub_tree = root.get_child("solver"); }
catch (std::exception &e)
    { std::cout << e.what() << std::endl; }

    analytic_corr = sub_tree.get<bool>("analytic corrections",false);
    theta = sub_tree.get<double>("theta",0.5);
    EPSILON = sub_tree.get<double>("epsilon",1e-16);
    DUMIN = sub_tree.get<double>("min(du)",1e-9);
    DUMAX = sub_tree.get<double>("max(du)",0.02);
    DTMIN = sub_tree.get<double>("min(dt)",1e-14);
    DTMAX = sub_tree.get<double>("max(dt)",1e-7);
    TAUR = 100.*DTMAX; // pourquoi 100 ?
    dt = sub_tree.get<double>("initial dt",1e-9);
    MAXITER = sub_tree.get<int>("max(iter)",500);
    REFRESH_PRC = sub_tree.get<int>("refresh prc every",20);
}
