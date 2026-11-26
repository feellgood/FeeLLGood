/** \file main.cpp
 * fg-merge-sol: combine a mesh file with a feeLLGood .sol file into a tab separated values
 * file giving the physical quantities as a function of position.
 */

#include <fstream>
#include <iostream>
#include <stdexcept>

#include <QCommandLineParser>
#include <QCoreApplication>

#include "gmshMesh.h"
#include "pyRepr.h"
#include "textTable.h"

namespace
    {
const char *description =
        "Combine a mesh file (which has the node coordinates) with a feeLLGood .sol file\n"
        "(which has the node's reduced magnetizations or electric potentials). Output the\n"
        ".sol physical quantities (magnetization or potential) as a function of position,\n"
        "in a tab-separated-values file: (x, y, z, physical_quantity...).\n"
        "Example:\n"
        "fg-merge-sol my_mesh.msh my_simul.sol -o my_results.tsv\n"
        "\n"
        "will generate an output file named my_results.tsv, without node indices, replaced\n"
        "by physical node coordinates, and associated physical quantities from my_simul.sol\n"
        "input file.";

/** Write to out the comments of the .sol file, with the "## columns:" line changed to name
 * the coordinates instead of the node index, then one line per node: its coordinates and the
 * physical quantities of the .sol file, without the node index. Throws std::runtime_error if
 * the .sol file has fewer data lines than the mesh has nodes. */
void merge(const fg::Mesh &mesh, const fg::TextTable &sol, std::ostream &out)
    {
    for (const std::string &line : sol.comments)
        {
        std::string value;
        if (fg::metadataValue(line, "## columns:", value))
            {  // remove the label of the index column
            const std::size_t tab = line.find('\t');
            const std::string rest = (tab == std::string::npos) ? "" : line.substr(tab + 1);
            out << "## columns: x\ty\tz\t" << rest << '\n';
            }
        else
            {
            out << line << '\n';
            }
        }
    if (sol.rows.size() < mesh.points.size())
        {
        throw std::runtime_error("the .sol file has fewer data lines than the mesh has nodes");
        }
    for (std::size_t i = 0; i < mesh.points.size(); ++i)
        {
        const fg::Vec3 &p = mesh.points[i];
        out << fg::pyRepr(p[0]) << '\t' << fg::pyRepr(p[1]) << '\t' << fg::pyRepr(p[2]);
        const std::vector<double> &row = sol.rows[i];
        for (std::size_t k = 1; k < row.size(); ++k)
            {  // column 0 is the node index
            out << '\t' << fg::pyRepr(row[k]);
            }
        out << '\n';
        }
    }
    }  // namespace

/** Parse the command line, read the mesh and the .sol file, and write the merged data to the
 * output file or to the standard output. Returns 1 on error. */
int main(int argc, char *argv[])
    {
    QCoreApplication app(argc, argv);
    QCoreApplication::setApplicationName("fg-merge-sol");
    QCoreApplication::setApplicationVersion(FG_VERSION);

    QCommandLineParser parser;
    parser.setApplicationDescription(description);
    parser.addHelpOption();
    parser.addVersionOption();
    parser.addOption({{"o", "output"}, "name of the output file, defaults to stdout", "output"});
    parser.addPositionalArgument("mesh_file", "name of the mesh file, in gmsh format");
    parser.addPositionalArgument("sol_file", "name of a matching feeLLGood .sol file");
    parser.process(app);
    const QStringList args = parser.positionalArguments();
    if (args.size() != 2)
        {
        std::cerr << "error: expected a mesh file and a .sol file\n";
        parser.showHelp(1);
        }

    try
        {
        fg::GmshSession gmshSession;
        const fg::Mesh mesh = fg::readMesh(args[0].toStdString(), true);
        const fg::TextTable sol = fg::readTextTable(args[1].toStdString());
        if (parser.isSet("output"))
            {
            const std::string outName = parser.value("output").toStdString();
            std::ofstream out(outName);
            if (!out)
                {
                throw std::runtime_error("cannot write to '" + outName + "'");
                }
            merge(mesh, sol, out);
            }
        else
            {
            merge(mesh, sol, std::cout);
            }
        }
    catch (const std::exception &e)
        {
        std::cerr << "error: " << e.what() << '\n';
        return 1;
        }
    return 0;
    }
