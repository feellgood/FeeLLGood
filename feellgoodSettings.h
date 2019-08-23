/** \file feellgoodSettings.h
\brief many settings to give some parameters to the solver, boundary conditions for the problem, the output file format wanted by the user. This is done mainly with the class Settings. 
*/

#include <iostream>
#include <vector>
#include <map>
#include <cmath>
#include <string>

#include "exprtk.hpp"

#include "pt3D.h"
#include "tetra.h"
#include "facette.h"

#ifndef feellgoodSettings_h
#define feellgoodSettings_h


/** \class Settings
container class to store many setting parameters, such as file names, parameters for the solver, output file format. It also handles text user interation through terminal, and some parsing functions. 
*/
class Settings{
	public:
	inline Settings() 
        {
        withTsv=true;
        withVtk = false;
        theta=0.5;
        analytic_corr = true;
        MAXITER = 500;
        REFRESH_PRC = 20;
        recentering = true;
        
        solverNbTh = 8;
        scalfmmNbTh = 8;
        
        s_table.add_variable("x",x);
        s_table.add_variable("y",y);
        s_table.add_variable("z",z);
        s_table.add_constants();
        } /**< default constructor */
	
	void infos(void);/**< some prints sent to terminal  */	
	
	void read(std::string fileJson);/**< read settings from a json file */
	
	inline void setPbName(std::string str) {pbName = str;} /**< setter for .msh file name  */
	inline std::string getPbName() {return pbName;}/**< getter for problem file name */

	inline void setSimName(std::string str) {simName = str;} /**< setter for .sol output file name  */
	inline std::string getSimName() {return simName;}/**< getter for output file name */

	inline void setScale(double s) {_scale = s;} /**< setter for geometrical scaling factor for physical coordinates of the mesh  */	
	inline double getScale(void) {return _scale;}/**< getter for geometrical scaling factor for physical coordinates of the mesh */ 

	/** maximum number of iterations setter for bicgstab */
    inline void set_MAXITER(int i) {MAXITER = i;} 
    
    /** refresh diagonal pre-conditionner every REFRESH_PRC for bicgstab */
    inline void set_REFRESH_PRC(int i) {REFRESH_PRC = i;}
	
	bool withTsv;/**< boolean flag to mention if you want output in txt tsv file format  */
	bool withVtk;/**< boolean flag to mention if you want output in txt vtk file format (readable by paraview) */
	double tf;/**< end time of the simulation */
	double dt;/**< step time of the simulation */

	int n1;/**< energy saved every n1 iterations */
	int n2;/**< magnetic configuration saved every n2 iterations */
	bool restore;/**< usefull to run a simulation from a previous calculation  */
    bool recentering;/**< to recenter DW or not */
	
	int solverNbTh;/**< nb of threads for the finite element solver */
	int scalfmmNbTh;/**< nb of threads for the computation of the demag field with scalfmm */
    
    triple Bext;/**< applied external field */
    
	std::string sMx;/**< string for analytical definition of Mx */
	std::string sMy;/**< string for analytical definition of My */
    std::string sMz;/**< string for analytical definition of Mz */
    
    std::string restoreFileName;/**< input file name for continuing a calculation (sol.in) */
	
	/** for time integration \f$ \theta \f$ scheme  */
	double theta;
	
    /** analytical corrections to the potential, extra contribution from surface charges */
    bool analytic_corr;
    
/** \f$ \epsilon \f$ is a small value, used to modify slightly J in (tet|facette).integrales */
    double EPSILON;

/** minimum value for du step */
    double DUMIN;	//1e-6 ; 1e-9 en dynamique

/** maximum value for du step */
    double DUMAX; // 0.1 en stat; 0.02 en dynamique

/** minimum step time for time integrator */
    double DTMIN; //1e-14;

/** maximum step time for time integrator */
    double DTMAX;//  1e-5 en stat ;  1e-7 en dynamique;

/** reduced \f$ \tau_r = 100 DT_{\mathrm{max}} \f$ */
    double TAUR;//    = 100.*DTMAX;
    
    int MAXITER;/**< maximum number of iteration for biconjugate gradient algorithm */
    
    int REFRESH_PRC;/**< refresh every REFRESH_PRC the diagonal preconditioner */
    
	/** this vector contains the material parameters for all regions for all the tetrahedrons */
	std::vector<Tetra::prm> paramTetra;

	/** \return index of the region in volume region container  */
	inline int findTetraRegionIdx(int r /**< [in] */) 
	{ 
	std::vector<Tetra::prm>::iterator result = std::find_if(paramTetra.begin(),paramTetra.end(),[r](Tetra::prm const& p){return(p.reg == r); }  ); 
	if (result == paramTetra.end()) return -1;
	else {return std::distance(paramTetra.begin(),result);}	
	};
	
	/** this vector contains the material parameters for all regions for all the facettes */
	std::vector<Facette::prm> paramFacette;

	/** relative path for output files (to be implemented) */
	std::string r_path_output_dir;
	
	/** \return index of the region in surface region container  */
	inline int findFacetteRegionIdx(int r /**< [in] */) 
	{ 
	std::vector<Facette::prm>::iterator result = std::find_if(paramFacette.begin(),paramFacette.end(),[r](Facette::prm const& p){return(p.reg == r); }  ); 
	if (result == paramFacette.end()) return -1;
	else {return std::distance(paramFacette.begin(),result);}	
	};
	
    /** parser magnetization compiler */
    inline void doCompile(void)
        {
        exprtk::parser<double> parser;
            
        expr_Mx.register_symbol_table(s_table);
        expr_My.register_symbol_table(s_table);
        expr_Mz.register_symbol_table(s_table);
        parser.compile(sMx,expr_Mx);
        parser.compile(sMy,expr_My);
        parser.compile(sMz,expr_Mz);
        };
    
    /** evaluation of the magnetization components through math expression */
    inline Pt::pt3D getValue(const Pt::pt3D &p)
        {
        x= p.x();
        y= p.y();
        z= p.z();
        
        double mx = expr_Mx.value();
        double my = expr_My.value();
        double mz = expr_Mz.value();
        Pt::pt3D mag = Pt::pt3D(mx,my,mz);
        mag.normalize();
        return mag;
        }
    
	private:
	
	double _scale;/**< scaling factor from gmsh files to feellgood */
	std::string simName;/**< simulation name */
	std::string pbName;     /**< mesh file, gmsh file format */
	
	double x;/**< working variable x for exprtk parser function f(x,y,z) */
	double y;/**< working variable y for exprtk parser function f(x,y,z) */
    double z;/**< working variable z for exprtk parser function f(x,y,z) */
	
	exprtk::symbol_table<double> s_table; /**< symbol table for exprtk math parser  */
	exprtk::expression<double> expr_Mx;/**< Mx expression */
    exprtk::expression<double> expr_My;/**< My expression */
    exprtk::expression<double> expr_Mz;/**< Mz expression */
};

#endif /* feellgoodSettings_h */
