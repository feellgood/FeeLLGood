#ifndef feellgoodSettings_h
#define feellgoodSettings_h

/** \file feellgoodSettings.h
\brief many settings to give some parameters to the solver, boundary conditions for the problem, the
output file format wanted by the user. This is done mainly with the class Settings.
*/

#include <cmath>
#include <iostream>
#include <map>
#include <string>
#include <vector>

#include <yaml-cpp/yaml.h>

#include "expression_parser.h"
#include "facette.h"
#include "spinTransferTorque.h"
#include "tetra.h"
#include "time_integration.h"

/** \class Settings
container class to store many setting parameters, such as file names, parameters for the solver,
output file format. It also handles text user interation through terminal, and some parsing
functions.
*/
class Settings
    {
public:
    /** default constructor */
    Settings();

    /** print out the YAML document defining the default settings */
    static void dumpDefaults();

    /** some prints sent to terminal */
    void infos(void);

    /** build a metadata string for .evol file */
    std::string evolMetadata(std::string realWorldTime) const;

    /** build a metadata string for .sol text files */
    std::string solMetadata(double t, std::string columnsTitle) const;

    /** read settings from a parsed YAML document */
    void read(YAML::Node);

    /** read settings from a YAML file. Returns true if a non-empty configuration is found. */
    bool read(std::string filename);

    /** returns numeric precision for .sol output text files */
    inline int getPrecision(void) const { return precision; }

    /** getter for fileDisplayName */
    inline std::string getFileDisplayName(void) const { return fileDisplayName; }

    /** setter for fileDisplayName */
    inline void setFileDisplayName(std::string _s) { fileDisplayName = _s; }

    /** setter for .msh file name */
    inline void setPbName(std::string str) { pbName = str; }

    /** getter for problem file name */
    inline std::string getPbName(void) const { return pbName; }

    /** setter for .sol output file name */
    inline void setSimName(std::string str) { simName = str; }

    /** getter for output file name */
    inline std::string getSimName(void) const { return simName; }

    /** setter for geometrical scaling factor for physical coordinates of the mesh */
    inline void setScale(const double s) { _scale = s; }

    /** getter for geometrical scaling factor for physical coordinates of the mesh */
    inline double getScale(void) const { return _scale; }

    /** maximum number of iterations setter for bicgstab */
    inline void set_MAXITER(int i) { MAXITER = i; }

    /** boolean flag to mention if you want output in txt tsv file format */
    bool withTsv;

    /** verbosity level, defaults to zero */
    int verbose;

    /** energy saved every time_step */
    double time_step;

    /** magnetic configuration saved every save_period time steps */
    int save_period;

    /** to recenter magnetization distribution or not */
    bool recenter;

    /** recentering direction, should be IDX_X|IDX_Y|IDX_Z */
    Nodes::index recentering_direction;

    /** threshold value to recenter or not versus avg(M_recentering_direction) */
    double threshold;

    /** nb of threads for the finite element solver */
    int solverNbTh;

    /** nb of threads for the computation of the demag field with scalfmm */
    int scalfmmNbTh;

    /** spin transfert torque parameters */
    STT p_stt;

    /** if spin transfer torque p_stt is fully initialized and boundary conditions ok, stt_flag is set to true */
    bool stt_flag;

    /** string for analytical definition of Mx */
    std::string sMx;

    /** string for analytical definition of My */
    std::string sMy;

    /** string for analytical definition of Mz */
    std::string sMz;

    /** string for a JavaScript function defining M, alternative to (sMx, sMy, sMz) */
    std::string sM;

    /** string for analytical definition of Bx */
    std::string sBx;

    /** string for analytical definition of By */
    std::string sBy;

    /** string for analytical definition of Bz */
    std::string sBz;

    /** string for a JavaScript function defining B, alternative to (sBx, sBy, sBz) */
    std::string sB;

    /** input file name for continuing a calculation (sol.in) */
    std::string restoreFileName;

    /** maximum value for du step */
    double DUMAX;  // 0.1 for magnetostatic simulations; 0.02 for the dynamics

    /** solver tolerance (eigen bicgstab) */
    double TOL;

    /** maximum number of iteration for biconjugate gradient algorithm */
    int MAXITER;

    /** ILU preconditioner dropping tolerance
    https://eigen.tuxfamily.org/dox/classEigen_1_1IncompleteLUT.html
    */
    double ILU_tol;

    /** ILU preconditioner filling factor
    https://eigen.tuxfamily.org/dox/classEigen_1_1IncompleteLUT.html
    */
    int ILU_fill_factor;

    /** this vector contains the material parameters for all regions for all the tetrahedrons */
    std::vector<Tetra::prm> paramTetra;

    /** \return index of the region in volume region container  */
    inline int findTetraRegionIdx(const std::string name /**< [in] */) const
        {
        std::vector<Tetra::prm>::const_iterator result =
                std::find_if(paramTetra.begin(), paramTetra.end(),
                             [name](Tetra::prm const &p) { return (p.regName == name); });

        int idx(-2);
        if (result == paramTetra.end())
            {
            idx = -1;
            }
        else
            {
            idx = std::distance(paramTetra.begin(), result);
            }
        return idx;
        };

    /** this vector contains the material parameters for all regions for all the facettes */
    std::vector<Facette::prm> paramFacette;

    /** relative path for output files (to be implemented) */
    std::string r_path_output_dir;

    /** contain the value names of the columns the user want in .evol file */
    std::vector<std::string> evol_columns;

    /** final integration time */
    double tf;

    /** minimal time step */
    double dt_min;

    /** maximal time step */
    double dt_max;

    /** \return index of the region in surface region container  */
    inline int findFacetteRegionIdx(const std::string name /**< [in] */) const
        {
        std::vector<Facette::prm>::const_iterator result =
                std::find_if(paramFacette.begin(), paramFacette.end(),
                             [name](Facette::prm const &p) { return (p.regName == name); });
        int idx(-2);

        if (result == paramFacette.end())
            {
            idx = -1;
            }
        else
            {
            idx = std::distance(paramFacette.begin(), result);
            }
        return idx;
        };

    /** evaluation of the magnetization components through math expression, each component of the
     * magnetization is a function of (x,y,z). It is not safe to call this method simultaneously
     * from multiple threads.
     * \return unit vector
     */
    inline Eigen::Vector3d getMagnetization(const Eigen::Ref<Eigen::Vector3d> p) const
        {
        Eigen::Vector3d tmp = mag_parser.get_vector(p);
        tmp.normalize();
        return tmp;
        }

    /** evaluation of the field components through math expression, each component of the field is a
     * function of (t).
     */
    inline Eigen::Vector3d getField(const double t_val) const
        {
        return nu0 * (field_parser.get_vector(t_val));
        }

private:
    int precision;               /**< numeric precision for .sol output text files */
    std::string fileDisplayName; /**< parameters file name : either a yaml file or standard input */
    double _scale;               /**< scaling factor from gmsh files to feellgood */
    std::string simName;         /**< simulation name */
    std::string pbName;          /**< mesh file, gmsh file format */
    VectorParser mag_parser;     /**< parser for the magnetization expressions */
    VectorParser field_parser;   /**< parser for the time dependant applied field expressions */
    };

#endif /* feellgoodSettings_h */
