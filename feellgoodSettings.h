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
#include "tetra.h"
#include "time_integration.h"

/** Specify how the initial magnetization is defined as a JavaScript function. It is either a
 * function of (x, y, z), the position coordinates of the node, or it takes a fourth parameter,
 * which is an array of strings with the names of the mesh volume regions the node belongs to.
 *
 */
enum mag_exprType
    {
    POSITION_ONLY,        ///< M(x, y, z)
    POSITION_AND_REGIONS  ///< M(x, y, z, regions)
    };

/** Specify how the applied magnetic field is defined using JavaScript functions. The choices are:
 *
 * * `RtoR3`: **B**(t), a vector function of time only, for a uniform applied field.
 * * `R4toR3`: **B**(x, y, z, t) = A(t) **f**(x, y, z), the product of a scalar function of space
 *   and a vector function of time.
 *
 * It is not possible to define the field as a general function **B**(x, y, z, t), as that would
 * require an unreasonable number of evaluations.
 */
enum field_exprType
    {
    UNDEF = -1,  ///< undefined
    RtoR3 = 1,   ///< **B**(t)
    R4toR3 = 2   ///< **B**(x, y, z, t) = A(t) **f**(x, y, z), where A(t) is a scalar
                 /// and **f**(x, y, z) is a vector
    };

/**
 * \class Settings
 *
 * \brief Container for all the settings provided by the user, with conversions to/from YAML.
 *
 * This class stores all the settings the user provides in their YAML configuration file. In
 * addition to storing the settings, the class can convert them to and from YAML text:
 *
 * - The two read() methods initialize the settings from the provided file or YAML::Node. Any
 *   setting not explicitly initialized keeps its default value.
 * - The toYaml() method serializes the settings as a YAML document which it prints to stdout.
 *
 * The default settings are provided in the source file default-settings.yml. During the build, the
 * contents of this file is embedded into the executable as a character array. This class'
 * constructor initializes the settings by parsing this embedded document.
 *
 * Javascript expressions (initial magnetization and applied field) are stored both in source form
 * (sM* and sB* members) and as private ExpressionParser instances. This class provides methods for
 * evaluating those expressions: getMagnetization(), getField(), getFieldSpace() and getFieldTime().
 * All these methods temporarily modify their corresponding ExpressionParser instances, it is thus
 * **not safe** to call them from multiple threads simultaneously.
 */
class Settings
    {
public:
    /** Constructor. The default values for data members are defined in the source file
     * default-settings.yml. A cmake rule invokes the linker `ld` in order to convert this file to
     * an object file which stores the original file's contents as a (not NUL-terminated) character
     * array. This constructor calls `read(YAML::Load(...))` on this embedded document in order to
     * initialize all data members.
     */
    Settings();

    /** Print to stdout the embedded YAML document defining the default settings. */
    static void dumpDefaults();

    /** Serialize these settings as a YAML document, and print them to stdout. */
    void toYaml(void);

    /** build a metadata string for .evol file */
    std::string evolMetadata() const;

    /** build a metadata string for .sol text files */
    std::string solMetadata(const double t) const;

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

    /** nb of threads for the computation of the demag field with scalfmm */
    int scalfmmNbTh;

    /** if true creates an output file of the solution of the electrostatic problem  */
    bool V_file;

    /** if spin_acc is true then spin accumulation contributions are computed */
    bool spin_acc;

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

    /** string for a JavaScript function defining B space: (x,y,z)->[,,] */
    std::string sB_space;

    /** string for a JavaScript function defining B time, (t) -> real */
    std::string sB_time;

    /** input file name for continuing a calculation (sol.in) */
    std::string restoreFileName;

    /** initial time for the simulation, if defined in the settings file, NAN otherwise */
    double initial_time = NAN;

    /** maximum value for du step */
    double DUMAX;  // 0.1 for magnetostatic simulations; 0.02 for the dynamics

    /** solver tolerance (eigen bicgstab) */
    double TOL;

    /** maximum number of iteration for biconjugate gradient algorithm */
    int MAXITER;

    /** this vector contains the material parameters for all regions for all the tetrahedrons */
    std::vector<Tetra::prm> paramTetra;

    /** \return index of the region in volume region container  */
    inline int findTetraRegionIdx(const std::string &name /**< [in] */) const
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

    /** Returns the initial magnetization type, which can be either `POSITION_ONLY` (depends only
     * on the position) or `POSITION_AND_REGIONS` (depends on the position and volume regions).
     */
    inline mag_exprType getMagType() const
        {
        if (mag_parser.parameter_count() < 4) return POSITION_ONLY;
        return POSITION_AND_REGIONS;
        }

    /** Evaluate the initial magnetization vector as a function of the position vector. This method
     * is **not thread safe**.
     * \return reduced magnetization (unit vector)
     */
    inline Eigen::Vector3d getMagnetization(const Eigen::Ref<Eigen::Vector3d> p) const
        {
        Eigen::Vector3d tmp = mag_parser.get_vector(p);
        tmp.normalize();
        return tmp;
        }

    /** Evaluate the initial magnetization vector as a function of the position vector and the set
     * of regions this node belongs to. This method is **not thread safe**.
     * \return reduced magnetization (unit vector)
     */
    inline Eigen::Vector3d getMagnetization(const Eigen::Ref<Eigen::Vector3d> p,
            const std::vector<std::string> &regions) const
        {
        Eigen::Vector3d tmp = mag_parser.get_vector(p, regions);
        tmp.normalize();
        return tmp;
        }

    /** Evaluate the applied field as a function of time. This is allowed only if field_type is
     * `RtoR3`. This method is **not thread safe**.
     * \return applied **H** field (not **B**!)
     */
    inline Eigen::Vector3d getField(const double t_val) const
        {
        return (field_parser.get_vector(t_val))/mu0;
        }

    /** Returns the field type, which can be either `RtoR3` (vector function of time) or `R4toR3`
     * (product of a time-dependent scalar and a position-dependent vector).
     */
    inline field_exprType getFieldType(void) const
        { return field_type; }

    /** Evaluate the “spatial” part of the applied field, which is a position-dependent vector. This
     * is allowed only if field_type is `R4toR3`. This method is **not thread safe**.
     * \return spatial part of applied **H** field (not **B**!)
     */
    inline Eigen::Vector3d getFieldSpace(const Eigen::Ref<Eigen::Vector3d> p) const
        { return field_space_parser.get_vector(p); }

    /** Evaluate the “time” part of the applied field, which is a time-dependent scalar. This is
     * allowed only if field_type is `R4toR3`. This method is **not thread safe**.
     * \return time part of applied **H** field (not **B**!)
     */
    inline double getFieldTime(const double t_val) const
        { return (field_time_parser.get_scalar(t_val))/mu0; }

private:
    using MetadataItem = std::pair<std::string, std::string>;  /**< type of userMetadata items */

    std::vector<MetadataItem> userMetadata;  /**< user-provided metadata for the output files */
    int precision;               /**< numeric precision for .sol output text files */
    std::string fileDisplayName; /**< parameters file name : either a yaml file or standard input */
    double _scale;               /**< scaling factor from gmsh files to feellgood */
    std::string simName;         /**< simulation name */
    std::string pbName;          /**< mesh file, gmsh file format */
    ExpressionParser mag_parser;     /**< parser for the magnetization expressions R³->R³ */

    /** applied field might be either expression defined (t) -> (Bx(t),By(t),Bz(t))
    or a couple of math expressions (t) -> A(t) and (x,y,z) -> (fx(x,y,z),fy(x,y,z),fz(x,y,z))
    with B(x,y,z,t) = A(t) * (fx(x,y,z),fy(x,y,z),fz(x,y,z))
    */
    field_exprType field_type;

    /** parser for the time dependant applied field expressions R->R³ */
    ExpressionParser field_parser;

    // field(x,y,z,t) = field_time(t) * field_space(x,y,z) with field_time R -> R and field_space R³->R³
    /** parser for the field time dependant scalar expression */
    ExpressionParser field_time_parser;

    /** parser for the field space dependant vector expressions */
    ExpressionParser field_space_parser;

    /** Metadata common to .evol and .sol files */
    std::ostringstream commonMetadata() const;

    /** to normalize input vector */
    const bool NORMALIZE = true;
    };

#endif /* feellgoodSettings_h */
