#include <iostream>
#include <unistd.h> // for getpid()
#include <signal.h>

#include "Utils/FTic.hpp"//for counter from scalfmm

#include "fem.h"
#include "mesh.h"
#include "linear_algebra.h"
#include "fmm_demag.h"
#include "time_integration.h"

#include "electrostatSolver.h"

// Catch SIGTERM in order to save the state before quitting.
volatile sig_atomic_t received_sigterm = 0;

static void sigterm_handler(int)
{
received_sigterm = 1;
}

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
struct Option
    {
    std::string short_opt, long_opt;
    bool *setting;
    };
struct Option options[] =
    {
    {"-v", "--verbose", &settings.verbose},
    {"", "", nullptr}  // sentinel
    };

int optind;
for (optind = 1; optind < argc; optind++)
    {
    char *opt = argv[optind];
    Option *o;
    for (o = options; o->setting; o++)
        {
        if (opt == o->short_opt || opt == o->long_opt)
            {
            *o->setting = true;
            break;
            }
        }
    if (!o->setting)  // option not found
        break;
    }
if (optind != argc - 1)
    {
    std::cerr << "Usage: feellgood [options] config_file.yml\n";
    exit(1);
    }
if (settings.verbose)
    std::cout << "Verbose mode active.\n";
std::string filename = argv[optind];
std::cout << "Loading settings from " << filename << "\n";
return filename;
}

int main(int argc,char* argv[])
{
Settings mySettings;
FTic counter;

prompt();

std::string filename = parseOptions(mySettings,argc,argv);
mySettings.read(filename);
timing t_prm = timing(mySettings.tf, mySettings.dt_min, mySettings.dt_max);
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
    electrostatSolver pot_solver = electrostatSolver(fem.msh,mySettings.p_stt,1e-8,mySettings.verbose,5000); 
    chronoElec.tac();
    std::cout << "\tSTT total computing time: " << chronoElec.elapsed() << " s\n" << std::endl;
    }
else
    { std::cout << "no spin transfer torque" << std::endl; }

scal_fmm::fmm myFMM(fem.msh,mySettings.verbose,mySettings.scalfmmNbTh);

// Catch SIGTERM.
struct sigaction action;
action.sa_handler = sigterm_handler;
action.sa_flags = 0;
sigemptyset(&action.sa_mask);
if (sigaction(SIGTERM, &action, NULL) == -1) {
    perror("SIGTERM");
    return EXIT_FAILURE;
}

int nt = time_integration(fem,mySettings,linAlg,myFMM,t_prm);
        
counter.tac();
std::cout << "\n  * iterations: " << nt << "\n  * total computing time: " << counter.elapsed() << " s\n" << std::endl;

return 0;
}
