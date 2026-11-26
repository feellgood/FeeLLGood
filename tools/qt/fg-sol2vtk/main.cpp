/** \file main.cpp
 * fg-sol2vtk: convert a feeLLGood .sol file and its gmsh mesh to a ParaView .vtu file.
 */

#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

#include <QCommandLineParser>
#include <QCoreApplication>
#include <QFileInfo>

#include <vtkFloatArray.h>
#include <vtkNew.h>
#include <vtkPointData.h>
#include <vtkPoints.h>
#include <vtkUnstructuredGrid.h>
#include <vtkXMLUnstructuredGridWriter.h>

#include "gmshMesh.h"
#include "textTable.h"

namespace
    {
const char *epilogue =
        "Build a .vtu file from feeLLGood output .sol file and .msh mesh file.\n"
        "By default the .vtu will contain magnetization components as a vector field.\n"
        "Example:\n"
        "  fg-sol2vtk myMesh.msh myMagConfig.sol\n"
        "If some columns are specified by their indices the .vtu file is built using them.\n"
        "Example:\n"
        "  fg-sol2vtk myMesh.msh myMagConfig.sol 1 2 3 4\n"
        "Would create a .vtu with four scalar fields, mx, my, mz and demagnetizing "
        "potential.";

/** Format a list of column indices as Python prints a list: [1, 2, 3]. */
std::string listToString(const std::vector<int> &indices)
    {
    std::string s = "[";
    for (std::size_t i = 0; i < indices.size(); ++i)
        {
        s += (i ? ", " : "") + std::to_string(indices[i]);
        }
    return s + "]";
    }

/** Build a float array named name (unnamed if empty) with one tuple per data line of sol,
 * made of the columns of the given indices. Throws std::runtime_error if an index is out of
 * the range of a line. */
vtkNew<vtkFloatArray> selectColumns(const fg::TextTable &sol, const std::vector<int> &indices,
                                    const std::string &name)
    {
    vtkNew<vtkFloatArray> array;
    array->SetNumberOfComponents(static_cast<int>(indices.size()));
    if (!name.empty())
        {
        array->SetName(name.c_str());
        }
    std::vector<double> tuple(indices.size());
    for (const auto &row : sol.rows)
        {
        for (std::size_t k = 0; k < indices.size(); ++k)
            {
            const int idx = indices[k];
            if (idx < 0 || static_cast<std::size_t>(idx) >= row.size())
                {
                throw std::runtime_error("index out of bounds: index = " + std::to_string(idx));
                }
            tuple[k] = row[idx];
            }
        array->InsertNextTuple(tuple.data());
        }
    return array;
    }
    }  // namespace

/** Parse the command line, build an unstructured grid of the tetrahedrons of the mesh with
 * the selected .sol columns as point data, and write it in binary .vtu format. Three columns
 * are written as a vector field, any other number as scalars. Returns 1 on error. */
int main(int argc, char *argv[])
    {
    QCoreApplication app(argc, argv);
    QCoreApplication::setApplicationName("fg-sol2vtk");
    QCoreApplication::setApplicationVersion(FG_VERSION);

    QCommandLineParser parser;
    parser.setApplicationDescription(
            QString("convert .sol to .vtu using gmsh mesh and .sol file\n\n") + epilogue);
    parser.addHelpOption();
    parser.addVersionOption();
    parser.addOption({"name", "name of the data in the output .vtu file", "name"});
    parser.addOption({{"o", "output"}, "output file name, extension .vtu will be added",
                      "output"});
    parser.addPositionalArgument("mshFileName", "a mesh file in gmsh format with .msh "
                                                "extension");
    parser.addPositionalArgument("solFileName", "a .sol feellgood file");
    parser.addPositionalArgument("N", "0-based indices of the columns to copy into output "
                                      ".vtu file", "[N...]");
    parser.process(app);
    const QStringList args = parser.positionalArguments();
    if (args.size() < 2)
        {
        std::cerr << "error: expected a mesh file and a .sol file\n";
        parser.showHelp(1);
        }

    const QString mshFileName = args[0];
    const QString solFileName = args[1];
    std::vector<int> indices;
    for (int i = 2; i < args.size(); ++i)
        {
        bool ok;
        indices.push_back(args[i].toInt(&ok));
        if (!ok)
            {
            std::cerr << "error: invalid column index '" << args[i].toStdString() << "'\n";
            return 1;
            }
        }
    if (indices.empty())
        {
        indices = {1, 2, 3};
        }

    const QString ext = "." + QFileInfo(mshFileName).suffix().toLower();
    if (ext != ".msh")
        {
        std::cerr << "extension " << ext.toStdString() << " not recognized.\n";
        return 1;
        }

    try
        {
        fg::GmshSession gmshSession;
        const fg::Mesh mesh = fg::readMesh(mshFileName.toStdString(), true);
        const fg::TextTable sol = fg::readTextTable(solFileName.toStdString());
        if (sol.rows.size() != mesh.points.size())
            {
            throw std::runtime_error("the .sol file has " + std::to_string(sol.rows.size())
                                     + " data lines, the mesh has "
                                     + std::to_string(mesh.points.size()) + " nodes");
            }
        vtkNew<vtkFloatArray> data =
                selectColumns(sol, indices, parser.value("name").toStdString());

        vtkNew<vtkPoints> points;
        points->SetNumberOfPoints(static_cast<vtkIdType>(mesh.points.size()));
        for (std::size_t i = 0; i < mesh.points.size(); ++i)
            {
            points->SetPoint(static_cast<vtkIdType>(i), mesh.points[i].data());
            }
        vtkNew<vtkUnstructuredGrid> ugrid;
        ugrid->SetPoints(points);
        ugrid->Allocate(static_cast<vtkIdType>(mesh.tetrahedrons.size()));
        for (const auto &t : mesh.tetrahedrons)
            {
            vtkIdType ids[4] = {static_cast<vtkIdType>(t[0]), static_cast<vtkIdType>(t[1]),
                                static_cast<vtkIdType>(t[2]), static_cast<vtkIdType>(t[3])};
            ugrid->InsertNextCell(VTK_TETRA, 4, ids);
            }

        std::string newFileName;
        if (parser.isSet("output"))
            {
            newFileName = parser.value("output").toStdString();
            }
        else
            {
            const std::string s = solFileName.toStdString();
            newFileName = s.substr(0, s.rfind('.'));
            }
        newFileName += ".vtu";

        // the values are on the nodes, so they are point data, not cell data
        if (indices.size() == 3)
            {
            std::cout << newFileName << " saved as 3D vectors, using " << listToString(indices)
                      << " columns\n";
            ugrid->GetPointData()->SetVectors(data);
            }
        else
            {
            std::cout << newFileName << " saved as " << indices.size() << "D scalars, using "
                      << listToString(indices) << " column(s)\n";
            ugrid->GetPointData()->SetScalars(data);
            }

        vtkNew<vtkXMLUnstructuredGridWriter> writer;
        writer->SetFileName(newFileName.c_str());
        writer->SetInputData(ugrid);
        writer->SetDataModeToBinary();
        if (writer->Write() != 1)
            {
            throw std::runtime_error("cannot write '" + newFileName + "'");
            }
        }
    catch (const std::exception &e)
        {
        std::cerr << "error: " << e.what() << '\n';
        return 1;
        }
    return 0;
    }
