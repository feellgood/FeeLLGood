/** \file feellgoodSettings.h
\brief many settings to give some parameters to the solver, boundary conditions for the problem, the output file format wanted by the user. This is done mainly with the class Settings. 
*/

#include <iostream>
#include <vector>
#include <map>
#include <cmath>
#include <string>

#ifndef feellgoodSettings_h
#define feellgoodSettings_h


typedef double triple[3];/**< a 3D point */

/** 
\return square of a number \f$ x^2 \f$
*/
inline double sq(double x /**< [in] */ ) {return x*x;}

/**
in place normalizing function of triple
a */
inline void normalize(triple &a /**< [in,out] */)
{
double norme=sqrt(sq(a[0])+sq(a[1])+sq(a[2]));
a[0]/= norme;a[1]/= norme;a[2]/= norme;
}

/** \struct Seq
Seq describe a sequence of field from \f$ B_{ini} \f$ to \f$ B_{fin} \f$ by steps \f$ dB \f$ with orientation \f$ \vec{a} \f$, should be a unit vecor.
*/    
struct Seq{
double Bini;/**< starting value */
double Bfin;/**< ending value */
double dB;/**< step field */
triple a;/**< direction (should be normalized) */
};

/* for statics 
const double EPSILON = 1e-40;
const double DUMAX   = 0.1;
const double DTMIN   = 1e-14;
const double DTMAX   = 1.e-5;
*/



/* for dynamics */

/** \f$ \epsilon \f$ is smallest possible value (for what value ? to check) */
const double EPSILON = 1e-40;

/** min for du step */
const double DUMIN   = 1e-9;	//1e-6

/** max for du step */
const double DUMAX   = 0.02;	// 0.02

/** minimum step time for time integrator */
const double DTMIN   = 1e-14;


/** maximum step time for time integrator */
const double DTMAX   = 1e-7;	// 1e-7

/** reduced /f$ /tau /f$ */
const double TAUR    = 100.*DTMAX;

/* physical constants */

const double mu0 = 4.*M_PI*1e-7;/**< \f$ \mu_0 = 4 \pi 10^{-7} \f$ */
const double nu0 = 1./mu0;/**< \f$ \nu_0 = 1/\mu_0 \f$ */

/** \class Settings
container class to store many setting parameters, such as file names, parameters for the solver, output file format. It also handles text user interation through terminal, and some parsing functions. 
*/
class Settings{
	public:
	inline Settings() { withTsv=true; withVtk = false;} /**< default constructor */
	
	void helloDialog();/**< kinda hello world! :-) */
	void printToTerminal(std::vector<Seq> &seq);/**< some prints sent to terminal  */	
	void dialog(std::vector<Seq> &seq);/**< interaction with user from terminal to fix some parameters of the simulation to run */
	

	inline void setPbName(std::string str) {pbName = str;} /**< setter for .msh file name  */
	inline std::string getPbName() {return pbName;}/**< getter for problem file name */

	inline void setSimName(std::string str) {simName = str;} /**< setter for .sol output file name  */
	inline std::string getSimName() {return simName;}/**< getter for output file name */

	inline void setScale(double s) {_scale = s;} /**< setter for geometrical scaling factor for physical coordinates of the mesh  */	
	inline double getScale(void) {return _scale;}/**< getter for geometrical scaling factor for physical coordinates of the mesh */ 

	bool withTsv;/**< boolean flag to mention if you want output in txt tsv file format  */
	bool withVtk;/**< boolean flag to mention if you want output in txt vtk file format (readable by paraview) */
	double tf;/**< end time of the simulation */
	double dt;/**< step time of the simulation */

	std::map < std::pair<std::string,int>,double>  param;/**< convenient map for association of parameters to keys, it stores different parameters of a simulation, such as suppression of the charges on some surfaces or not, the material parameter values as \f$ Ms \f$ or \f$ \alpha \f$ */
	int n1;/**< energy saved every n1 iterations */
	int n2;/**< magnetic configuration saved every n2 iterations */
	int restore;/**< usefull to run a simulation from a previous calculation  */

	private:
	
	void extract_comment(std::istream &flux);/**< parser */

	double _scale;/**< scaling factor from gmsh files to feellgood */
	std::string simName;/**< simulation name */
	std::string pbName;     /**< mesh file, gmsh file format */
};

#endif /* feellgoodSettings_h */
