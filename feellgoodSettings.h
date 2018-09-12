/** \file feellgoodSettings.h
\brief many settings to give some parameters to the solver, boundary conditions for the problem, the output file format wanted by the user. This is done mainly with the class Settings. 
*/

#include <iostream>
#include <vector>
#include <map>
#include <cmath>
#include <string>

#include "tetra.h"
#include "facette.h"

#ifndef feellgoodSettings_h
#define feellgoodSettings_h



/** \struct Seq
Seq describe a sequence of field from \f$ B_{ini} \f$ to \f$ B_{fin} \f$ by steps \f$ dB \f$ with orientation \f$ \vec{a} \f$, should be a unit vector.
*/    
struct Seq{
double Bini;/**< starting value */
double Bfin;/**< ending value */
double dB;/**< step field */
triple a;/**< field direction componants */
};



/** \class Settings
container class to store many setting parameters, such as file names, parameters for the solver, output file format. It also handles text user interation through terminal, and some parsing functions. 
*/
class Settings{
	public:
	inline Settings() { withTsv=true; withVtk = false; theta=0.5;} /**< default constructor */
	
	void helloDialog();/**< kinda hello world! */
	void printToTerminal(std::vector<Seq> &seq);/**< some prints sent to terminal  */	
	void dialog(std::vector<Seq> &seq);/**< text interaction with user from terminal to fix some parameters of the simulation to run */
	
	void read(std::vector<Seq> &seq);/**< read settings from a .json file */
	
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

	std::map < std::pair<std::string,int>,double>  param;/**< map for association of parameters to keys, it stores different parameters of a simulation, such as suppression of the charges on some surfaces or not, the material parameter values as \f$ Ms \f$ or \f$ \alpha \f$ */
	int n1;/**< energy saved every n1 iterations */
	int n2;/**< magnetic configuration saved every n2 iterations */
	bool restore;/**< usefull to run a simulation from a previous calculation  */

	/** for time integration \f$ \theta \f$ scheme  */
	double theta;
	
	/** this vector contains the material parameters for all regions for all the tetrahedrons */
	std::vector<Tetra::prm> paramTetra;

	/** this vector contains the material parameters for all regions for all the facettes */
	std::vector<Facette::prm> paramFacette;

	private:
	
	void extract_comment(std::istream &flux);/**< parser */

	double _scale;/**< scaling factor from gmsh files to feellgood */
	std::string simName;/**< simulation name */
	std::string pbName;     /**< mesh file, gmsh file format */
};

#endif /* feellgoodSettings_h */
