#include <iostream>
#include <unistd.h> // for getpid()

#include "Utils/FTic.hpp"//for counter from scalfmm

#include "linear_algebra.h"
#include "fmm_demag.h"

using namespace std;

void prompt(void)
{
std::cout << "\n\t ┌────────────────────────────┐\n";
std::cout <<   "\t │         FeeLLGood          │\n";
std::cout <<   "\t │        version " << feellgood_version << "       │\n";
std::cout <<   "\t │      cnrs Grenoble-INP     │\n";
std::cout <<   "\t └────────────────────────────┘\n";
std::cout << "\t process\t\t" << getpid() << std::endl;
}

int time_integration(Fem &fem,Settings &settings /**< [in] */,LinAlgebra &linAlg /**< [in] */,scal_fmm::fmm &myFMM  /**< [in] */);

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
	
mySettings.read(fileJson);
mySettings.infos();

Fem fem = Fem(mySettings);

counter.tic();

fem.infos();

//once fem containers are ok, linAlgebra object is built
LinAlgebra linAlg(mySettings,fem.node,fem.tet,fem.fac);
linAlg.set_Hext(fem.Hext[0],fem.Hext[1],fem.Hext[2]);

scal_fmm::fmm myFMM = scal_fmm::fmm(fem,mySettings.verbose,mySettings.scalfmmNbTh);

myFMM.calc_demag(fem,mySettings);
   
int nt = time_integration(fem,mySettings,linAlg,myFMM);
        
counter.tac();
cout << "\n  * iterations: " << nt << "\n  * total computing time: " << counter.elapsed() << " s\n--- the end ---\n" << endl;
    
return 0;
}
