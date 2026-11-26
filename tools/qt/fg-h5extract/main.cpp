/** \file main.cpp
 * fg-h5extract: extract from a HDF5 file written by fg-sol2hdf5 the mesh, the settings, the
 * .evol file or the .sol files of a simulation.
 */

#include <cstdio>
#include <fstream>
#include <iostream>
#include <stdexcept>

#include <QCommandLineParser>
#include <QCoreApplication>
#include <QFileInfo>

#include "h5Utils.h"

namespace
    {
const std::string EMPTY_GROUP = "emptyGroupName";
const std::string EMPTY_BASE = "emptyBaseName";

/** Write to fileName the text metadata followed by the text content. Throws
 * std::runtime_error if the file cannot be written. */
void writeText(const std::string &fileName, const std::string &content,
               const std::string &metadata)
    {
    std::ofstream f(fileName, std::ios::binary);
    if (!f)
        {
        throw std::runtime_error("cannot write to '" + fileName + "'");
        }
    f << metadata << content;
    }

/** Write to fileName the metadata, an empty line, then the rows of the table with tab
 * separated values in %.18e format, as numpy.savetxt does. Throws std::runtime_error if the
 * file cannot be written. */
void writeTable(const std::string &fileName, const std::vector<std::vector<double>> &rows,
                const std::string &metadata)
    {
    std::ofstream f(fileName, std::ios::binary);
    if (!f)
        {
        throw std::runtime_error("cannot write to '" + fileName + "'");
        }
    f << metadata << '\n';
    char buf[32];
    for (const auto &row : rows)
        {
        for (std::size_t k = 0; k < row.size(); ++k)
            {
            std::snprintf(buf, sizeof(buf), "%.18e", row[k]);
            f << (k ? "\t" : "") << buf;
            }
        f << '\n';
        }
    }

/** Rebuild the comment lines of a .sol file from the attributes of its dataset. */
std::string solMetadata(const H5::DataSet &ds)
    {
    std::string s = fg::readStringAttribute(ds, "comments");
    s += "## feeLLGood version: " + fg::readStringAttribute(ds, "feeLLGood version") + '\n';
    s += "## settings file: " + fg::readStringAttribute(ds, "settings file") + '\n';
    s += "## time: " + fg::attributeToString(ds.openAttribute("time")) + '\n';
    s += "## columns: " + fg::readStringAttribute(ds, "columns") + '\n';
    return s;
    }

/** Write the .sol file of the dataset named name of group to fileName. */
void writeSolFile(const H5::Group &group, const std::string &name, const std::string &fileName)
    {
    const H5::DataSet ds = group.openDataSet(name);
    writeTable(fileName, fg::readTable(ds), solMetadata(ds));
    }

/** Write the .evol file of group to baseName.evol. */
void extractEvol(const H5::Group &group, const std::string &baseName)
    {
    const H5::DataSet ds = group.openDataSet("evolFile");
    writeTable(baseName + ".evol", fg::readTable(ds), fg::readStringAttribute(ds, "metadata"));
    }

/** Write the mesh of group to baseName.msh. */
void extractMesh(const H5::Group &group, const std::string &baseName)
    {
    writeText(baseName + ".msh", fg::readStringAttribute(group, "mesh"), "");
    }

/** Write the settings of group to baseName.yml. */
void extractSettings(const H5::Group &group, const std::string &baseName)
    {
    writeText(baseName + ".yml", fg::readStringAttribute(group, "settings"), "");
    }

/** Write every .sol dataset of group (all but evolFile) to baseName<dataset name>.sol. */
void extractDatasets(const H5::Group &group, const std::string &baseName)
    {
    for (const std::string &name : fg::memberNames(group))
        {
        if (name != "evolFile")
            {
            writeSolFile(group, name, baseName + name + ".sol");
            }
        }
    }

/** Print the names of the attributes of a group, then its datasets with their attributes
 * and their values. */
void viewGroup(const H5::H5File &file, const std::string &groupName)
    {
    std::cout << groupName << '\n';
    const H5::Group group = file.openGroup(groupName);
    for (const std::string &a : fg::attributeNames(group))
        {
        std::cout << "  a: " << a << '\n';
        }
    for (const std::string &d : fg::memberNames(group))
        {
        std::cout << "  d: " << d << '\n';
        const H5::DataSet ds = group.openDataSet(d);
        for (const std::string &a : fg::attributeNames(ds))
            {
            std::cout << "    a: " << a << " : " << fg::attributeToString(ds.openAttribute(a))
                      << '\n';
            }
        }
    }

/** Print the content of the group named groupName, or of all the groups if it is not
 * given. */
void view(const H5::H5File &file, const std::string &groupName)
    {
    if (groupName != EMPTY_GROUP && file.nameExists(groupName))
        {
        viewGroup(file, groupName);
        }
    else
        {
        for (const std::string &g : fg::memberNames(file))
            {
            viewGroup(file, g);
            }
        }
    }
    }  // namespace

/** Parse the command line, then view the content of the file or extract from a group what
 * --what asks for: evol, mesh, settings, datasets, the name of one dataset, or everything by
 * default. Returns 1 on error. */
int main(int argc, char *argv[])
    {
    QCoreApplication app(argc, argv);
    QCoreApplication::setApplicationName("fg-h5extract");
    QCoreApplication::setApplicationVersion(FG_VERSION);

    QCommandLineParser parser;
    parser.setApplicationDescription(
            "extract from hdf5 mesh, feellgood settings(yaml), .evol or .sol file(s)\n\n"
            "Extracts data from hdf5 file and writes file(s) from hdf5 feellgood input file.");
    parser.addHelpOption();
    parser.addVersionOption();
    parser.addOption({{"g", "groupname"}, "group name where to extract", "groupname",
                      QString::fromStdString(EMPTY_GROUP)});
    parser.addOption({{"b", "basename"}, "a base name for output file(s)", "basename",
                      QString::fromStdString(EMPTY_BASE)});
    parser.addOption({"view", "view in standard output the group content"});
    parser.addOption({{"w", "what"},
                      "what to extract: evol, mesh, settings, datasets or a dataset name; "
                      "everything if not given",
                      "what", "undefined"});
    parser.addPositionalArgument("h5FileName", "an input hdf5 file name");
    parser.process(app);
    const QStringList args = parser.positionalArguments();
    if (args.size() != 1)
        {
        std::cerr << "error: expected one hdf5 file name\n";
        parser.showHelp(1);
        }

    const std::string fileName = args[0].toStdString();
    const std::string groupName = parser.value("groupname").toStdString();
    const std::string baseName = parser.value("basename").toStdString();
    const std::string what = parser.value("what").toStdString();

    if (!QFileInfo::exists(args[0]))
        {
        std::cerr << "error: file " << fileName << " not found\n";
        return 1;
        }

    try
        {
        H5::Exception::dontPrint();
        const H5::H5File file(fileName, H5F_ACC_RDONLY);
        const bool groupExists = groupName != EMPTY_GROUP && file.nameExists(groupName);

        if (parser.isSet("view") && (groupName == EMPTY_GROUP || groupExists))
            {
            view(file, groupName);
            return 0;
            }
        if (!groupExists)
            {
            std::cerr << "error: group " << groupName << " does not exist in file " << fileName
                      << '\n';
            return 1;
            }

        const H5::Group group = file.openGroup(groupName);
        if (what == "evol")
            {
            extractEvol(group, baseName);
            }
        else if (what == "mesh")
            {
            extractMesh(group, baseName);
            }
        else if (what == "settings")
            {
            extractSettings(group, baseName);
            }
        else if (what == "datasets")
            {
            extractDatasets(group, baseName);
            }
        else if (group.nameExists(what))
            {
            writeSolFile(group, what, baseName + what + ".sol");
            }
        else
            {  // extract all
            extractEvol(group, baseName);
            extractMesh(group, baseName);
            extractSettings(group, baseName);
            extractDatasets(group, baseName);
            }
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
