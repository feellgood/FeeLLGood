#include <iostream>
#include <unistd.h> // for getpid()

#include "Utils/FTic.hpp"//for counter from scalfmm

#include "linear_algebra.h"
#include "fmm_demag.h"

#include "time_integration.h"

using namespace std;


int time_integration(Fem &fem,Settings &settings /**< [in] */,LinAlgebra &linAlg /**< [in] */,scal_fmm::fmm &myFMM  /**< [in] */,timing &t_prm);

void prompt(void)
{
std::cout << "\n\t ┌────────────────────────────┐\n";
std::cout <<   "\t │         FeeLLGood          │\n";
std::cout <<   "\t │        version " << feellgood_version << "       │\n";
std::cout <<   "\t │      cnrs Grenoble-INP     │\n";
std::cout <<   "\t └────────────────────────────┘\n";
std::cout << "\t process\t\t" << getpid() << std::endl;
}

int main(int argc,char* argv[])
{
Settings mySettings;

FTic counter;
 
string fileJson;

if(argc<2)
	{
	fileJson = "settings.json";
	cout << "using default settings from " << fileJson << " JSON file." <<endl; 
	}
else 
	{
	fileJson = argv[1]; // argv[0] is "./feellgood"
	cout << "using loaded settings from " << fileJson << " JSON file." <<endl;	
	}

prompt();

timing t_prm;
mySettings.read(t_prm,fileJson);
Fem fem = Fem(mySettings,t_prm);

mySettings.infos();
t_prm.infos();
fem.infos();

counter.tic();
//once fem containers are ok, linAlgebra object is built
LinAlgebra linAlg(mySettings,fem.node,fem.tet,fem.fac);

scal_fmm::fmm myFMM = scal_fmm::fmm(fem,mySettings.verbose,mySettings.scalfmmNbTh);

myFMM.calc_demag(fem,mySettings);
std::cout << "first calc demag done." << std::endl;

int nt = time_integration(fem,mySettings,linAlg,myFMM,t_prm);
        
counter.tac();
cout << "\n  * iterations: " << nt << "\n  * total computing time: " << counter.elapsed() << " s\n--- the end ---\n" << endl;
    
return 0;
}
