#include <iostream>
#include <unistd.h> // for getpid()

#include "Utils/FTic.hpp"//for counter from scalfmm

#include "fem.h"
#include "mesh.h"
#include "linear_algebra.h"
#include "fmm_demag.h"
#include "time_integration.h"

#include "electrostatSolver.h"

int time_integration(Fem &fem,Settings &settings /**< [in] */,LinAlgebra &linAlg /**< [in] */,scal_fmm::fmm &myFMM  /**< [in] */,timing &t_prm);

inline std::string spaceString(int nbSpace)
{std::string S;
    for(int i=0;i<nbSpace;i++) {S += " ";}
return S;
}

void prompt(void)
{
    std::stringstream SsId;
    SsId << getpid();
std::cout << "\n\t ┌────────────────────────────┐\n";
std::cout <<   "\t │         FeeLLGood          │\n";
std::cout <<   "\t │        version " << feellgood_version << spaceString(28-16-feellgood_version.length() ) <<"│\n";
std::cout <<   "\t │      cnrs Grenoble-INP     │\n";
std::cout <<   "\t │    feellgood.neel.cnrs.fr  │\n";
std::cout <<   "\t ├────────────────────────────┤\n";
std::cout <<   "\t │      process "<< SsId.str() <<  spaceString(28-14-SsId.str().length() ) <<"│\n";
std::cout <<   "\t └────────────────────────────┘\n";
}

std::string parseOptions(Settings &settings,int argc,char* argv[])
{
std::string fileJson;
std::cout << "parsing options..." << std::endl;

if(argc<2)
	{
	std::cout << "no JSON file provided, see FeeLLGood online documentation to create some settings using Python script settingsMaker.py" << std::endl; 
	exit(1);
    }
else if (argc == 2)
        {
        fileJson = argv[1]; // argv[0] is "./feellgood"
        std::cout << "using loaded settings from " << fileJson << " JSON file.\n";	
        }
    else if ((argc == 3)&&( strcmp(argv[1],"-v") == 0))
        {
            std::cout << "verbose mode active\n"  << std::endl;
            settings.verbose = true;
            fileJson = argv[2];
        }
return fileJson;
}

int main(int argc,char* argv[])
{
Settings mySettings;
FTic counter;

prompt();

std::string fileJson = parseOptions(mySettings,argc,argv);
timing t_prm = mySettings.read(fileJson);
Fem fem = Fem(mySettings,t_prm);

if(mySettings.verbose)
    {
    mySettings.infos();
    t_prm.infos();
    fem.infos();
    }
    
std::cout<<"start computing..."<<std::endl;
counter.tic();
LinAlgebra linAlg(mySettings,fem.msh);

if (mySettings.stt_flag)
    {
    FTic chronoElec;
    chronoElec.tic();    
    electrostatSolver pot_solver = electrostatSolver(fem.msh,mySettings.p_stt,1e-8,mySettings.verbose,false,5000); 
    chronoElec.tac();
    std::cout << "\tSTT total computing time: " << chronoElec.elapsed() << " s\n" << std::endl;
    }
else
    { std::cout << "no spin transfer torque" << std::endl; }

scal_fmm::fmm myFMM(fem,mySettings.verbose,mySettings.scalfmmNbTh);

int nt = time_integration(fem,mySettings,linAlg,myFMM,t_prm);
        
counter.tac();
std::cout << "\n  * iterations: " << nt << "\n  * total computing time: " << counter.elapsed() << " s\n" << std::endl;

return 0;
}
