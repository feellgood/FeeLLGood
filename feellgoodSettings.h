#ifndef feellgoodSettings_h
#define feellgoodSettings_h

/** \file feellgoodSettings.h
\brief many settings to give some parameters to the solver, boundary conditions for the problem, the output file format wanted by the user. This is done mainly with the class Settings. 
*/

#include <iostream>
#include <vector>
#include <map>
#include <cmath>
#include <string>

#include <yaml-cpp/yaml.h>

#include "mag_parser.h"
#include "time_integration.h"
#include "tetra.h"
#include "facette.h"
#include "spinTransferTorque.h"

/** \class Settings
container class to store many setting parameters, such as file names, parameters for the solver, output file format. It also handles text user interation through terminal, and some parsing functions. 
*/
class Settings{
	public:
    /** default constructor */
    Settings();

    /** print out the YAML document defining the default settigs */
    static void dumpDefaults();
    
    /** some prints sent to terminal */	
    void infos(void);

    /** read settings from a parsed YAML document */
	void read(YAML::Node);
	
    /** read settings from a YAML file. Returns true if a non-empty configuration is found. */
	bool read(std::string filename);
	
	inline void setPbName(std::string str) {pbName = str;} /**< setter for .msh file name  */
	inline std::string getPbName(void) const {return pbName;}/**< getter for problem file name */

	inline void setSimName(std::string str) {simName = str;} /**< setter for .sol output file name  */
	inline std::string getSimName(void) const {return simName;}/**< getter for output file name */

	inline void setScale(const double s) {_scale = s;} /**< setter for geometrical scaling factor for physical coordinates of the mesh  */	
	inline double getScale(void) const {return _scale;}/**< getter for geometrical scaling factor for physical coordinates of the mesh */ 

	/** maximum number of iterations setter for bicgstab */
    inline void set_MAXITER(int i) {MAXITER = i;} 
    
    /** refresh diagonal pre-conditionner every REFRESH_PRC for bicgstab */
    inline void set_REFRESH_PRC(int i) {REFRESH_PRC = i;}
	
	bool withTsv;/**< boolean flag to mention if you want output in txt tsv file format  */
	bool withVtk;/**< boolean flag to mention if you want output in txt vtk file format (readable by paraview) */
	
	bool verbose;/**< if true the user gets many printing feedback in terminal */
	
	double time_step;/**< energy saved every time_step */
	int save_period;/**< magnetic configuration saved every save_period time steps */
	
    bool recenter;/**< to recenter magnetization distribution or not */
	Pt::index recentering_direction;/**< recentering direction, should be IDX_X|IDX_Y|IDX_Z */
	double threshold;/**< threshold value to recenter or not versus avg(M_recentering_direction) */
	
	int solverNbTh;/**< nb of threads for the finite element solver */
	int scalfmmNbTh;/**< nb of threads for the computation of the demag field with scalfmm */
    
    STT p_stt;/**< spin transfert torque parameters */
    
    bool stt_flag;/**< if spin transfer torque p_stt is fully initialized and boundary conditions ok, stt_flag is set to true */
    
    std::string sMx;/**< string for analytical definition of Mx */
	std::string sMy;/**< string for analytical definition of My */
    std::string sMz;/**< string for analytical definition of Mz */
    
    std::string sBx;/**< string for analytical definition of Bx */
	std::string sBy;/**< string for analytical definition of By */
    std::string sBz;/**< string for analytical definition of Bz */
    
    std::string restoreFileName;/**< input file name for continuing a calculation (sol.in) */

/** maximum value for du step */
    double DUMAX; // 0.1 en stat; 0.02 en dynamique
    
    int MAXITER;/**< maximum number of iteration for biconjugate gradient algorithm */
    
    int REFRESH_PRC;/**< refresh every REFRESH_PRC the diagonal preconditioner */
    
	/** this vector contains the material parameters for all regions for all the tetrahedrons */
	std::vector<Tetra::prm> paramTetra;

	/** \return index of the region in volume region container  */
	inline int findTetraRegionIdx(const std::string name /**< [in] */) const
	{
	std::vector<Tetra::prm>::const_iterator result = std::find_if(paramTetra.begin(),paramTetra.end(),[name](Tetra::prm const& p){return(p.regName == name); }  ); 
	
	int idx(-2);
	if (result == paramTetra.end()) 
		{idx = -1;}
	else
		{idx = std::distance(paramTetra.begin(),result);}	
	return idx;
	};
	
	/** this vector contains the material parameters for all regions for all the facettes */
	std::vector<Facette::prm> paramFacette;

	/** relative path for output files (to be implemented) */
	std::string r_path_output_dir;
	
    /** contain the value names of the columns the user want in .evol file */
    std::vector<std::string> evol_columns;
    
    /**  if true the first line of .evol file is the title of each column tsv format, starting with \# */
    bool evol_header;

    /* timing parameters */
    double tf;      /**< final integration time */
    double dt_min;  /**< minimal time step */
    double dt_max;  /**< maximal time step */
    
	/** \return index of the region in surface region container  */
	inline int findFacetteRegionIdx(const std::string name /**< [in] */) const
	{ 
	std::vector<Facette::prm>::const_iterator result = std::find_if(paramFacette.begin(),paramFacette.end(),[name](Facette::prm const& p){return(p.regName == name); }  ); 
	int idx(-2);
	
	if (result == paramFacette.end())
		{ idx = -1; }
	else
		{ idx = std::distance(paramFacette.begin(),result);}	
	return idx;
	};
	
    /** parser magnetization compiler */
    inline void doCompile3Dprm(void)
        { mag_parser.set_expressions(sMx, sMy, sMz); };
    
    /** parser time dependant field compiler */
    inline void doCompile1Dprm(void)
        { field_parser.set_expressions(sBx, sBy, sBz); };
    
        
    /** evaluation of the magnetization components through math expression, each component of the magnetization is a function of (x,y,z). 
     \return unit vector
     */
    inline Pt::pt3D getValue(const Pt::pt3D &p) const
        { return mag_parser.get_magnetization(p); }
    
    /** evaluation of the field components through math expression, each component of the field is a function of (t). 
     */
    inline Pt::pt3D getValue(const double t_val) const
        { return nu0*(field_parser.get_timeDepField(t_val)); }
    
    
	private:
	
	double _scale;/**< scaling factor from gmsh files to feellgood */
	std::string simName;/**< simulation name */
	std::string pbName;     /**< mesh file, gmsh file format */
    MagnetizationParser mag_parser;  /**< parser for the magnetization expressions */
    TimeDepFieldParser field_parser; /**< parser for the time dependant applied field expressions */
};

#endif /* feellgoodSettings_h */
