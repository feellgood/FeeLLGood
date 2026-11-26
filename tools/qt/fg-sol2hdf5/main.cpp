/** \file main.cpp
 * fg-sol2hdf5: archive in a HDF5 file the .sol files of a simulation, its .evol file, its
 * settings file and its mesh.
 */

#include <fstream>
#include <iostream>
#include <map>
#include <sstream>
#include <stdexcept>

#include <QCommandLineParser>
#include <QCoreApplication>
#include <QDir>
#include <QFileInfo>

#include <yaml-cpp/yaml.h>

#include "h5Utils.h"
#include "textTable.h"

namespace
    {
const char *epilogue =
        "Build a hdf5 file from feellgood output .sol files and settings file.\n"
        "If no .sol file name is provided, all .sol files in the directory are considered.";

/** Metadata of a .sol file, stored as attributes of its dataset. */
struct SolMetadata
    {
    std::string comments;  //!< comment lines not listed below, with their end of line
    std::string version;   //!< value of "## feeLLGood version:"
    std::string settings;  //!< value of "## settings file:"
    double time = 0;       //!< value of "## time:"
    std::string columns;   //!< value of "## columns:"
    };

/** Read a whole text file. Throws std::runtime_error if it cannot be opened. */
std::string readFile(const std::string &fileName)
    {
    std::ifstream f(fileName, std::ios::binary);
    if (!f)
        {
        throw std::runtime_error("cannot open file '" + fileName + "'");
        }
    std::ostringstream ss;
    ss << f.rdbuf();
    return ss.str();
    }

/** Names of the files of the current directory matching pattern, in alphabetical order. */
QStringList filesMatching(const QString &pattern)
    {
    return QDir::current().entryList({pattern}, QDir::Files, QDir::Name);
    }

/** Name of the .evol file of the current directory: the first one with a warning if there
 * are several. Throws std::runtime_error if there is none. */
std::string inputEvolFile()
    {
    const QStringList files = filesMatching("*.evol");
    if (files.isEmpty())
        {
        throw std::runtime_error("no .evol file in the current directory");
        }
    if (files.size() > 1)
        {
        std::cout << "Warning: more than one .evol file, " << files[0].toStdString()
                  << " will be used\n";
        }
    return files[0].toStdString();
    }

/** Extract the metadata from the comment lines of a .sol file. */
SolMetadata solMetadata(const fg::TextTable &sol)
    {
    SolMetadata md;
    for (const std::string &line : sol.comments)
        {
        std::string value;
        if (fg::metadataValue(line, "## feeLLGood version:", value))
            {
            md.version += value;
            }
        else if (fg::metadataValue(line, "## settings file:", value))
            {
            md.settings += value;
            }
        else if (fg::metadataValue(line, "## time:", value))
            {
            md.time = std::stod(value);
            }
        else if (fg::metadataValue(line, "## columns:", value))
            {
            md.columns += value;
            }
        else
            {
            md.comments += line + "\n";
            }
        }
    return md;
    }

/** Name of the mesh file given by mesh.filename in the settings file, without its directory.
 * Throws if the settings file cannot be parsed or has no such entry. */
std::string meshFileNameFromSettings(const std::string &settingsFileName)
    {
    const YAML::Node settings = YAML::LoadFile(settingsFileName);
    const std::string meshFileName = settings["mesh"]["filename"].as<std::string>();
    const std::size_t idx = meshFileName.rfind('/');
    return (idx == std::string::npos) ? meshFileName : meshFileName.substr(idx + 1);
    }

/** Name of the dataset of a .sol file: "solFile" followed by the part of the file name from
 * its first '_' to ".sol", for instance solFile_iter100 for sim_iter100.sol. */
std::string datasetName(const std::string &solFileName)
    {
    const std::size_t begin = solFileName.find('_');
    const std::size_t end = solFileName.find(".sol");
    if (begin == std::string::npos || end == std::string::npos || end < begin)
        {
        return "solFile";
        }
    return "solFile" + solFileName.substr(begin, end - begin);
    }

/** Keep the first n columns of each row. Throws std::runtime_error if a row is shorter. */
std::vector<std::vector<double>> firstColumns(const std::vector<std::vector<double>> &rows,
                                              std::size_t first, std::size_t n,
                                              const std::string &fileName)
    {
    std::vector<std::vector<double>> out;
    out.reserve(rows.size());
    for (const auto &row : rows)
        {
        if (row.size() < first + n)
            {
            throw std::runtime_error("'" + fileName + "': too few columns");
            }
        out.emplace_back(row.begin() + first, row.begin() + first + n);
        }
    return out;
    }
    }  // namespace

/** Parse the command line and write the group of the HDF5 file: one dataset per .sol file
 * (columns 1 to 4, with the metadata as attributes), the evolFile dataset with the comments
 * of the .evol file as attribute, and the settings and the mesh as text attributes of the
 * group. Returns 1 on error. */
int main(int argc, char *argv[])
    {
    QCoreApplication app(argc, argv);
    QCoreApplication::setApplicationName("fg-sol2hdf5");
    QCoreApplication::setApplicationVersion(FG_VERSION);

    QCommandLineParser parser;
    parser.setApplicationDescription(
            QString("convert .sol to hdf5 using .sol file(s)\n\n") + epilogue);
    parser.addHelpOption();
    parser.addVersionOption();
    parser.addOption({{"a", "append"}, "append to an existing hdf5 file a group", "groupName"});
    parser.addPositionalArgument("outputFileName", "an output file name with h5 extension");
    parser.addPositionalArgument("settingsFileName", "a settings file");
    parser.addPositionalArgument("fileName", ".sol feellgood files", "[fileName...]");
    parser.process(app);
    const QStringList args = parser.positionalArguments();
    if (args.size() < 2)
        {
        std::cerr << "error: expected an output file and a settings file\n";
        parser.showHelp(1);
        }

    const std::string outputFileName = args[0].toStdString();
    const std::string settingsFileName = args[1].toStdString();
    QStringList solFileNames = args.mid(2);

    try
        {
        H5::Exception::dontPrint();
        const std::string meshFileName = meshFileNameFromSettings(settingsFileName);
        if (solFileNames.isEmpty())
            {
            solFileNames = filesMatching("*.sol");
            }
        // read everything before creating the output file, so that an error leaves no
        // partial file
        const std::string evolFileName = inputEvolFile();
        const fg::TextTable evol = fg::readTextTable(evolFileName);
        const std::string settingsText = readFile(settingsFileName);
        const std::string meshText = readFile(meshFileName);

        const bool append = parser.isSet("append");
        const std::string groupName = append ? parser.value("append").toStdString() : "results";
        const bool exists = QFileInfo::exists(QString::fromStdString(outputFileName));
        if (append && !exists)
            {
            std::cout << "New file " << outputFileName << " created, append to group "
                      << groupName << '\n';
            }
        if (!append && exists)
            {
            std::cout << "Warning: " << outputFileName << " already exists, overwriting.\n";
            }
        H5::H5File file(outputFileName,
                        append ? (exists ? H5F_ACC_RDWR : H5F_ACC_EXCL) : H5F_ACC_TRUNC);
        if (append && file.nameExists(groupName))
            {
            std::cerr << "Error: group name " << groupName << " already exists.\n";
            return 1;
            }
        H5::Group group = file.createGroup(groupName);

        for (const QString &solFileName : solFileNames)
            {
            const std::string name = solFileName.toStdString();
            const fg::TextTable sol = fg::readTextTable(name);
            const SolMetadata md = solMetadata(sol);
            H5::DataSet ds = fg::writeTable(group, datasetName(name),
                                            firstColumns(sol.rows, 1, 4, name));
            fg::writeStringAttribute(ds, "comments", md.comments);
            fg::writeStringAttribute(ds, "feeLLGood version", md.version);
            fg::writeStringAttribute(ds, "settings file", md.settings);
            fg::writeDoubleAttribute(ds, "time", md.time);
            fg::writeStringAttribute(ds, "columns", md.columns);
            }

        std::string evolMetadata;
        for (const std::string &line : evol.comments)
            {
            evolMetadata += line + "\n";
            }
        H5::DataSet evolDs = fg::writeTable(group, "evolFile", evol.rows);
        fg::writeStringAttribute(evolDs, "metadata", evolMetadata);

        fg::writeStringAttribute(group, "settings", settingsText);
        // the mesh is copied as a text blob
        fg::writeStringAttribute(group, "mesh", meshText);
        }
    catch (const H5::Exception &e)
        {
        std::cerr << "error: HDF5: " << e.getDetailMsg() << '\n';
        return 1;
        }
    catch (const std::exception &e)
        {
        std::cerr << "error: " << e.what() << '\n';
        return 1;
        }
    return 0;
    }
