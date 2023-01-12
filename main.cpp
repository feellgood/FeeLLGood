#include <iostream>
#include <unistd.h> // for getpid(), stat()
#include <stdio.h>  // for perror()
#include <sys/stat.h>  // for mkdir(), stat()
#include <sys/types.h>  // for mkdir(), stat()
#include <signal.h>

#include "Utils/FTic.hpp"//for counter from scalfmm

#include "fem.h"
#include "mesh.h"
#include "linear_algebra.h"
#include "fmm_demag.h"
#include "time_integration.h"

#include "electrostatSolver.h"

// Catch some deadly signals in order to save the state before quitting.
volatile sig_atomic_t received_signal = 0;

static void signal_handler(int signal_number)
{
received_signal = signal_number;
}

// Create the output directory if it does not exist yet.
static void create_dir_if_needed(std::string dirname)
{
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

int time_integration(Fem &fem,Settings &settings /**< [in] */,LinAlgebra &linAlg /**< [in] */,scal_fmm::fmm &myFMM  /**< [in] */,timing &t_prm);

// Return the number of characters in an UTF-8-encoded string.
static int char_length(const std::string &s)
{
    int len = 0;
    for (const char *p = s.c_str(); *p; ++p) {
        if ((*p & 0xc0) != 0x80)  // count bytes not matching 0x10xxxxxx
            len++;
    }
    return len;
}

// Return a string padded with spaces to a given length.
static std::string pad(const std::string &s, int length)
{
    length -= char_length(s);
    if (length < 0)
        length = 0;
    return s + std::string(length, ' ');
}

void prompt(void)
{
std::cout << "\t┌───────────────────────────────┐\n";
std::cout << "\t│           feeLLGood           │\n";
std::cout << "\t│      CNRS Grenoble – INP      │\n";
std::cout << "\t│ http://feellgood.neel.cnrs.fr │\n";
std::cout << "\t└───────────────────────────────┘\n";
}

std::string parseOptions(Settings &settings,
      int argc, char* argv[], unsigned int &random_seed)
{
bool print_help = false;
bool print_version = false;
bool print_defaults = false;
bool verify = false;
bool use_fixed_seed = false;
struct Option
    {
    std::string short_opt, long_opt;
    const char *help;
    bool *setting;
    };
struct Option options[] =
    {
    {"-h", "--help", "display short help and exit", &print_help},
    {"-V", "--version", "display version information and exit", &print_version},
    {"", "--print-defaults", "print default settings and exit", &print_defaults},
    {"", "--verify", "verify a settings file and exit", &verify},
    {"-v", "--verbose", "enable verbose mode", &settings.verbose},
    {"", "--seed", "set random seed", &use_fixed_seed},
    {"", "", nullptr, nullptr}  // sentinel
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
            if (o->long_opt == "--seed" && optind < argc - 1)
               {
               random_seed = atol(argv[++optind]);
               }
            break;
            }
        }
    if (!o->setting)  // option not found
        break;
    }
if (print_help)
    {
    std::cout << "Usage: " << argv[0] << " [options] settings_file\n";
    std::cout << "Options:\n";
    for (Option *o = options; o->setting; o++) {
        std::cout << "  ";
        if (!o->short_opt.empty())
            std::cout << o->short_opt << " ";
        else
            std::cout << "   ";
        std::cout << pad(o->long_opt, 17) << " " << o->help << "\n";
    }
    exit(0);
    }
if (print_version)
    {
    std::cout << "feeLLGood " << feellgood_version << "\n";
    exit(0);
    }
if (print_defaults)
    {
    Settings::dumpDefaults();
    exit(0);
    }
if (optind != argc - 1)
    {
    std::cerr << "Usage: feellgood [options] config_file.yml\n";
    exit(1);
    }
std::string filename = argv[optind];
if (verify)
    {
    settings.read(filename);
    settings.infos();
    exit(0);
    }
if (!use_fixed_seed)
   {
   // Use a truly random seed.
   std::random_device rd;
   random_seed = rd();
   }
return filename;
}

int main(int argc,char* argv[])
{
Settings mySettings;
FTic counter;

unsigned int random_seed;
std::string filename = parseOptions(mySettings, argc, argv, random_seed);
srand(random_seed);
prompt();
if (mySettings.verbose)
    std::cout << "verbose mode:      on\n";
std::cout << "feeLLGood version: " << feellgood_version << '\n';
std::cout << "process ID:        " << std::to_string(getpid()) << '\n';
std::cout << "random seed:       " << random_seed << '\n';
std::string fileDisplayName = filename=="-" ? "standard input" : filename;
std::cout << "settings file:     " << fileDisplayName << '\n';
if (!mySettings.read(filename))
    {
    std::cerr << "Error: no settings found.\n";
    return 1;
    }
std::cout << "mesh file:         " << mySettings.getPbName() << '\n';
std::cout << "output directory:  " << mySettings.r_path_output_dir << " ";
create_dir_if_needed(mySettings.r_path_output_dir);
timing t_prm = timing(mySettings.tf, mySettings.dt_min, mySettings.dt_max);
Fem fem = Fem(mySettings,t_prm);

if(mySettings.verbose)
    {
    std::cout << "-- settings: -----------------------------------\n";
    mySettings.infos();
    std::cout << "-- end of settings -----------------------------\n";
    t_prm.infos();
    fem.msh.infos();
    }
    
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

scal_fmm::fmm myFMM(fem.msh,mySettings.verbose,mySettings.scalfmmNbTh);

// Catch SIGINT and SIGTERM.
struct sigaction action;
action.sa_handler = signal_handler;
action.sa_flags = 0;
sigemptyset(&action.sa_mask);
if (sigaction(SIGINT, &action, NULL) == -1) {
    perror("SIGINT");
    return EXIT_FAILURE;
}
if (sigaction(SIGTERM, &action, NULL) == -1) {
    perror("SIGTERM");
    return EXIT_FAILURE;
}

int nt = time_integration(fem,mySettings,linAlg,myFMM,t_prm);
        
counter.tac();
double total_time = counter.elapsed();
std::cout << "\nComputing time:\n\n";
std::cout << "    total: " << total_time << " s";
if (total_time >= 3600)
    std::cout << " (" << total_time/3600 << " h)";
else if (total_time >= 60)
    std::cout << " (" << total_time/60 << " min)";
std::cout << '\n';
std::cout << "    per time step: " << total_time/nt << " s\n";
return 0;
}
