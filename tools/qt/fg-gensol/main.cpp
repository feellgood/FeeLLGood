/** \file main.cpp
 * fg-gensol: generate an initial magnetization file sol.in for feeLLGood from a mesh file and
 * a magnetization function per volume region, defined in JavaScript in init_cfg.js.
 */

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdio>
#include <ctime>
#include <fstream>
#include <iostream>
#include <limits>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

#include <QCommandLineParser>
#include <QCoreApplication>
#include <QFile>
#include <QFileInfo>
#include <QJSEngine>
#include <QJSValue>
#include <QSysInfo>
#include <QTextStream>

#include <gmsh.h>

#include "gmshMesh.h"
#include "jsMath.h"

namespace
    {
const char *epilogue =
        "fg-gensol generates initial magnetization configurations (sol.in files) for\n"
        "feeLLGood simulations from mesh files, in gmsh format. It detects the named\n"
        "physical volume regions of the mesh. On the first run, it writes init_cfg.js,\n"
        "where a JavaScript function of the node coordinates x, y, z returning\n"
        "[mx, my, mz] is defined for each region. Edit it, then run fg-gensol again to\n"
        "write sol.in. The functions of the Math object can be used without the Math.\n"
        "prefix, for instance: function (x, y, z) { return [cos(x), sin(x), 0]; }";

const QString configFile = "init_cfg.js";

/** A named physical volume region of the mesh. */
struct Region
    {
    std::string name;
    std::vector<std::size_t> nodes;  //!< indices of the nodes of its elements, sorted
    };

/** Named physical volume regions of the mesh currently loaded in gmsh, in the order of
 * their physical tags. Unnamed physical groups are ignored. Throws std::runtime_error if
 * there is no named volume region. */
std::vector<Region> volumeRegions(std::size_t nbNodes)
    {
    std::vector<std::pair<int, int>> groups;
    gmsh::model::getPhysicalGroups(groups, 3);
    std::vector<Region> regions;
    for (const auto &[dim, tag] : groups)
        {
        Region r;
        gmsh::model::getPhysicalName(dim, tag, r.name);
        if (r.name.empty())
            {
            continue;
            }
        std::vector<int> entities;
        gmsh::model::getEntitiesForPhysicalGroup(dim, tag, entities);
        for (int e : entities)
            {
            std::vector<int> types;
            std::vector<std::vector<std::size_t>> elemTags, nodeTags;
            gmsh::model::mesh::getElements(types, elemTags, nodeTags, dim, e);
            for (const auto &tags : nodeTags)
                {
                for (std::size_t t : tags)
                    {
                    r.nodes.push_back(fg::nodeIndex(t, nbNodes));
                    }
                }
            }
        std::sort(r.nodes.begin(), r.nodes.end());
        r.nodes.erase(std::unique(r.nodes.begin(), r.nodes.end()), r.nodes.end());
        regions.push_back(std::move(r));
        }
    if (regions.empty())
        {
        throw std::runtime_error("no named volume regions (dimension 3) found in mesh");
        }
    return regions;
    }

/** s as a single quoted JavaScript string literal. */
std::string jsString(const std::string &s)
    {
    std::string out = "'";
    for (char c : s)
        {
        if (c == '\'' || c == '\\')
            {
            out += '\\';
            }
        out += c;
        }
    return out + "'";
    }

/** JavaScript line defining the default magnetization [0, 0, 1] of a region. */
std::string defaultFunctionLine(const std::string &regionName)
    {
    return "    " + jsString(regionName) + ": function (x, y, z) { return [0, 0, 1]; },";
    }

/** Write init_cfg.js, defining the default magnetization [0, 0, 1] for every region, sorted
 * by name. Throws std::runtime_error if the file cannot be written. */
void autogenerateConfig(const std::vector<Region> &regions)
    {
    std::vector<std::string> names;
    for (const Region &r : regions)
        {
        names.push_back(r.name);
        }
    std::sort(names.begin(), names.end());

    std::ofstream f(configFile.toStdString());
    if (!f)
        {
        throw std::runtime_error("cannot write " + configFile.toStdString());
        }
    f << "// Auto-generated configuration file of fg-gensol\n"
         "// Edit the functions below to define the magnetization vector field of each "
         "region.\n"
         "// Each function takes the node coordinates x, y, z and returns [mx, my, mz].\n"
         "// The functions of the Math object can be used without the Math. prefix.\n"
         "\n"
         "({\n";
    for (const std::string &n : names)
        {
        f << defaultFunctionLine(n) << '\n';
        }
    f << "})\n";

    std::cout << "✅ File '" << configFile.toStdString() << "' created with " << regions.size()
              << " regions.\n"
              << "➡️  Edit " << configFile.toStdString()
              << " to define magnetization vector fields for each region.\n"
              << "➡️  Then run fg-gensol again to generate the sol.in file.\n";
    }

/** Evaluate init_cfg.js and return the object mapping region names to functions. Throws
 * std::runtime_error if the file cannot be read, has an error, or does not evaluate to an
 * object. */
QJSValue loadConfig(QJSEngine &engine)
    {
    QFile f(configFile);
    if (!f.open(QIODevice::ReadOnly | QIODevice::Text))
        {
        throw std::runtime_error("cannot read " + configFile.toStdString());
        }
    QTextStream in(&f);
    in.setCodec("UTF-8");
    const QJSValue config = engine.evaluate(in.readAll(), configFile);
    if (config.isError())
        {
        throw std::runtime_error(fg::jsErrorMessage(config).toStdString());
        }
    if (!config.isObject())
        {
        throw std::runtime_error("it must evaluate to an object mapping region names to "
                                 "functions");
        }
    return config;
    }

/** Magnetization of every node, computed with the function of its region and normalized.
 * Nodes shared by several regions get the value of the last one. Exits with a message if a
 * function fails, returns something else than three numbers, or if a node belongs to no
 * region. */
std::vector<fg::Vec3> applyRegionFunctions(const fg::Mesh &mesh,
                                           const std::vector<Region> &regions,
                                           const QJSValue &config)
    {
    const double nan = std::numeric_limits<double>::quiet_NaN();
    std::vector<fg::Vec3> m(mesh.points.size(), {nan, nan, nan});
    for (const Region &r : regions)
        {
        QJSValue func = config.property(QString::fromStdString(r.name));
        for (std::size_t node : r.nodes)
            {
            const fg::Vec3 &p = mesh.points[node];
            const QJSValue res = func.call({p[0], p[1], p[2]});
            if (res.isError())
                {
                throw std::runtime_error("region '" + r.name + "': "
                                         + fg::jsErrorMessage(res).toStdString());
                }
            if (!res.isArray() || res.property("length").toInt() != 3)
                {
                throw std::runtime_error("region '" + r.name
                                         + "': the function must return [mx, my, mz]");
                }
            for (quint32 k = 0; k < 3; ++k)
                {
                const QJSValue v = res.property(k);
                if (!v.isNumber())
                    {
                    throw std::runtime_error("region '" + r.name
                                             + "': the function must return numbers");
                    }
                m[node][k] = v.toNumber();
                }
            }
        }

    bool zero = false, uninitialized = false;
    for (fg::Vec3 &v : m)
        {
        const double n = std::sqrt(v[0] * v[0] + v[1] * v[1] + v[2] * v[2]);
        if (std::isnan(n))
            {
            uninitialized = true;
            }
        else if (n > 0)
            {
            v = {v[0] / n, v[1] / n, v[2] / n};
            }
        else
            {
            zero = true;
            }
        }
    if (zero)
        {
        std::cerr << "⚠️  m contains zero-norm vectors\n";
        }
    if (uninitialized)
        {
        std::cerr << "m contains non-initialized values\n";
        throw std::runtime_error("❌ abort");
        }
    return m;
    }

/** Local time in ISO 8601 format with microseconds, as Python's datetime.now().isoformat(). */
std::string isoNow()
    {
    using namespace std::chrono;
    const auto now = system_clock::now();
    const std::time_t t = system_clock::to_time_t(now);
    const long us = static_cast<long>(
            duration_cast<microseconds>(now.time_since_epoch()).count() % 1000000);
    std::tm tm{};
    localtime_r(&t, &tm);
    char buf[64];
    std::strftime(buf, sizeof(buf), "%Y-%m-%dT%H:%M:%S", &tm);
    std::string s = buf;
    if (us != 0)
        {
        std::snprintf(buf, sizeof(buf), ".%06ld", us);
        s += buf;
        }
    return s;
    }

/** Write the magnetization to fileName in the .sol format, with a zero potential column.
 * Throws std::runtime_error if the file cannot be written. */
void writeSolutionFile(const std::vector<fg::Vec3> &m, const std::string &fileName = "sol.in")
    {
    const std::vector<std::string> header = {
            "## generator: fg-gensol",
            "## hostname: " + QSysInfo::machineHostName().toStdString(),
            "## real-world time: " + isoNow(),
            "## time: 0",
            "## columns: idx\tmx\tmy\tmz\tphi"};
    std::ofstream f(fileName);
    if (!f)
        {
        throw std::runtime_error("cannot write " + fileName);
        }
    for (const std::string &line : header)
        {
        f << line << '\n';
        }
    char buf[128];
    for (std::size_t i = 0; i < m.size(); ++i)
        {
        std::snprintf(buf, sizeof(buf), "%zu\t%+.7e\t%+.7e\t%+.7e\t%.7e\n", i, m[i][0],
                      m[i][1], m[i][2], 0.0);
        f << buf;
        }
    std::cout << "✅ File '" << fileName << "' successfully written (" << header.size()
              << " header lines and " << m.size() << " data lines).\n";
    }
    }  // namespace

/** Parse the command line and read the mesh. Write init_cfg.js if it does not exist,
 * otherwise check that it defines every region and write sol.in. Returns 1 on error. */
int main(int argc, char *argv[])
    {
    QCoreApplication app(argc, argv);
    QCoreApplication::setApplicationName("fg-gensol");
    QCoreApplication::setApplicationVersion(FG_VERSION);

    QCommandLineParser parser;
    parser.setApplicationDescription(QString("generate a sol.in file\n\n") + epilogue);
    parser.addHelpOption();
    parser.addVersionOption();
    parser.addPositionalArgument("mshFileName", "a mesh file in gmsh format");
    parser.process(app);
    const QStringList args = parser.positionalArguments();
    if (args.size() != 1)
        {
        std::cerr << "error: expected one mesh file\n";
        parser.showHelp(1);
        }

    try
        {
        fg::GmshSession gmshSession;
        const fg::Mesh mesh = fg::readMesh(args[0].toStdString());
        const std::vector<Region> regions = volumeRegions(mesh.points.size());

        if (!QFileInfo::exists(configFile))
            {
            std::cout << "⚠️  File '" << configFile.toStdString() << "' not found.\n";
            autogenerateConfig(regions);
            return 0;
            }

        QJSEngine engine;
        fg::installMathGlobals(engine);
        QJSValue config;
        try
            {
            config = loadConfig(engine);
            }
        catch (const std::exception &e)
            {
            std::cerr << "⚠️  Error in '" << configFile.toStdString() << "': " << e.what()
                      << '\n';
            return 1;
            }

        std::vector<std::string> missing;
        for (const Region &r : regions)
            {
            if (!config.property(QString::fromStdString(r.name)).isCallable())
                {
                missing.push_back(r.name);
                }
            }
        if (!missing.empty())
            {
            std::cerr << "⚠️  Regions without a function in " << configFile.toStdString()
                      << ". Add these lines to it, then edit them:\n";
            for (const std::string &n : missing)
                {
                std::cerr << defaultFunctionLine(n) << '\n';
                }
            return 1;
            }

        std::cout << "✅ All regions are defined, generating sol.in...\n";
        writeSolutionFile(applyRegionFunctions(mesh, regions, config));
        }
    catch (const std::exception &e)
        {
        std::cerr << "error: " << e.what() << '\n';
        return 1;
        }
    return 0;
    }
