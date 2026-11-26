/** \file main.cpp
 * fg-plot: basic quiver rendering of a feeLLGood .sol file, in 3D or in a 2D subspace.
 */

#include <array>
#include <iostream>
#include <stdexcept>
#include <vector>

#include <QApplication>
#include <QCommandLineParser>
#include <QFileInfo>
#include <QJSEngine>
#include <QSurfaceFormat>

#include <QVTKOpenGLNativeWidget.h>
#include <vtkActor.h>
#include <vtkArrowSource.h>
#include <vtkCubeAxesActor.h>
#include <vtkDoubleArray.h>
#include <vtkGenericOpenGLRenderWindow.h>
#include <vtkGlyph3D.h>
#include <vtkNew.h>
#include <vtkPointData.h>
#include <vtkPoints.h>
#include <vtkPolyData.h>
#include <vtkPolyDataMapper.h>
#include <vtkProperty.h>
#include <vtkRenderer.h>
#include <vtkTextActor.h>
#include <vtkTextProperty.h>

#include "gmshMesh.h"
#include "jsMath.h"
#include "quiver2DWidget.h"
#include "textTable.h"

namespace
    {
/** Position and magnetization of a node: x, y, z, mx, my, mz. */
using Node = std::array<double, 6>;

/** Keep the nodes satisfying the predicate, a JavaScript expression of x, y, z, mx, my, mz
 * (Python's and, or, not are accepted), then one node out of int(100 / percent). An empty
 * predicate keeps all the nodes. Throws std::runtime_error if the predicate has an error or
 * no node satisfies it. */
std::vector<Node> filtering(const std::vector<Node> &nodes, const QString &predicate,
                            double percent)
    {
    std::vector<Node> kept;
    if (predicate.isEmpty())
        {
        kept = nodes;
        }
    else
        {
        QJSEngine engine;
        fg::installMathGlobals(engine);
        QJSValue pred = engine.evaluate("(function (x, y, z, mx, my, mz) { return ("
                                        + fg::pythonToJs(predicate) + "); })");
        if (pred.isError())
            {
            throw std::runtime_error("invalid filter: " + fg::jsErrorMessage(pred).toStdString());
            }
        for (const Node &n : nodes)
            {
            const QJSValue r = pred.call({n[0], n[1], n[2], n[3], n[4], n[5]});
            if (r.isError())
                {
                throw std::runtime_error("filter: " + fg::jsErrorMessage(r).toStdString());
                }
            if (r.toBool())
                {
                kept.push_back(n);
                }
            }
        if (kept.empty())
            {
            throw std::runtime_error("no node satisfies the filter " + predicate.toStdString());
            }
        }
    const std::size_t stride = static_cast<std::size_t>(100.0 / percent);
    std::vector<Node> out;
    for (std::size_t i = 0; i < kept.size(); i += stride)
        {
        out.push_back(kept[i]);
        }
    return out;
    }

/** Window with the 2D quiver plot of the components idx[2], idx[3] of the nodes at the
 * positions given by the components idx[0], idx[1]. */
QWidget *plot2D(const std::vector<Node> &nodes, const std::array<int, 4> &idx)
    {
    std::array<std::vector<double>, 4> c;
    for (const Node &n : nodes)
        {
        for (int k = 0; k < 4; ++k)
            {
            c[k].push_back(n[idx[k]]);
            }
        }
    return new Quiver2DWidget(c[0], c[1], c[2], c[3], "Magnetization subspace");
    }

/** Window with the 3D rendering of the magnetization of the nodes: blue arrows of length
 * qLength along the normalized magnetization, centered on the nodes, with x, y, z axes. */
QWidget *plot3D(const std::vector<Node> &nodes, double qLength)
    {
    vtkNew<vtkPoints> points;
    vtkNew<vtkDoubleArray> dirs;
    dirs->SetNumberOfComponents(3);
    for (const Node &n : nodes)
        {
        const double norm = std::sqrt(n[3] * n[3] + n[4] * n[4] + n[5] * n[5]);
        if (!(norm > 0))
            {
            continue;
            }
        const double d[3] = {n[3] / norm, n[4] / norm, n[5] / norm};
        points->InsertNextPoint(n[0] - 0.5 * qLength * d[0], n[1] - 0.5 * qLength * d[1],
                                n[2] - 0.5 * qLength * d[2]);
        dirs->InsertNextTuple(d);
        }
    vtkNew<vtkPolyData> poly;
    poly->SetPoints(points);
    poly->GetPointData()->SetVectors(dirs);

    vtkNew<vtkArrowSource> arrow;
    vtkNew<vtkGlyph3D> glyph;
    glyph->SetInputData(poly);
    glyph->SetSourceConnection(arrow->GetOutputPort());
    glyph->SetVectorModeToUseVector();
    glyph->OrientOn();
    glyph->SetScaleModeToDataScalingOff();
    glyph->SetScaleFactor(qLength);

    vtkNew<vtkPolyDataMapper> mapper;
    mapper->SetInputConnection(glyph->GetOutputPort());
    mapper->ScalarVisibilityOff();
    vtkNew<vtkActor> actor;
    actor->SetMapper(mapper);
    actor->GetProperty()->SetColor(0.122, 0.467, 0.706);  // matplotlib's default blue

    vtkNew<vtkRenderer> renderer;
    renderer->SetBackground(1, 1, 1);
    renderer->AddActor(actor);

    vtkNew<vtkCubeAxesActor> axes;
    axes->SetBounds(poly->GetBounds());
    axes->SetCamera(renderer->GetActiveCamera());
    axes->SetXTitle("x");
    axes->SetYTitle("y");
    axes->SetZTitle("z");
    axes->SetFlyModeToOuterEdges();
    axes->SetGridLineLocation(vtkCubeAxesActor::VTK_GRID_LINES_FURTHEST);
    axes->DrawXGridlinesOn();
    axes->DrawYGridlinesOn();
    axes->DrawZGridlinesOn();
    for (int k = 0; k < 3; ++k)
        {
        axes->GetTitleTextProperty(k)->SetColor(0, 0, 0);
        axes->GetLabelTextProperty(k)->SetColor(0, 0, 0);
        }
    axes->GetXAxesLinesProperty()->SetColor(0, 0, 0);
    axes->GetYAxesLinesProperty()->SetColor(0, 0, 0);
    axes->GetZAxesLinesProperty()->SetColor(0, 0, 0);
    renderer->AddActor(axes);

    vtkNew<vtkTextActor> title;
    title->SetInput("Magnetization");
    title->GetTextProperty()->SetColor(0, 0, 0);
    title->GetTextProperty()->SetFontSize(18);
    title->GetTextProperty()->SetJustificationToCentered();
    title->GetTextProperty()->SetVerticalJustificationToTop();
    title->GetPositionCoordinate()->SetCoordinateSystemToNormalizedViewport();
    title->SetPosition(0.5, 0.97);
    renderer->AddActor2D(title);
    renderer->ResetCamera();

    auto *widget = new QVTKOpenGLNativeWidget;
    vtkNew<vtkGenericOpenGLRenderWindow> window;
    widget->setRenderWindow(window);
    window->AddRenderer(renderer);
    widget->resize(800, 800);
    return widget;
    }
    }  // namespace

/** Parse the command line, read the mesh and the .sol file, filter the nodes and show the 2D
 * or the 3D plot. Returns 1 on error. */
int main(int argc, char *argv[])
    {
    QSurfaceFormat::setDefaultFormat(QVTKOpenGLNativeWidget::defaultFormat());
    QApplication app(argc, argv);
    QApplication::setApplicationName("fg-plot");
    QApplication::setApplicationVersion(FG_VERSION);

    QCommandLineParser parser;
    parser.setApplicationDescription(
            "make a basic quiver rendering of a .sol file\n\n"
            "Render a .sol file with a .msh and a corresponding .sol file.");
    parser.addHelpOption();
    parser.addVersionOption();
    parser.addOption({{"P", "percent"}, "Percentage of the field vectors rendered", "percent",
                      "100"});
    parser.addOption({{"q", "qlength"}, "quiver length", "qlength", "1.0"});
    parser.addOption({{"f", "filter"},
                      "predicate: an expression with x,y,z,mx,my,mz returning a boolean, in "
                      "JavaScript syntax; and, or, not are accepted, Math functions can be "
                      "used without the Math. prefix",
                      "filter"});
    parser.addOption({{"s", "subspace"},
                      "2D quiver plot of x,y,z,mx,my,mz read from input files, user has to "
                      "provide four indices, example \"0 1 3 4\" to plot on positions (x,y) "
                      "the set of quivers (mx,my)",
                      "subspace"});
    parser.addPositionalArgument("meshFileName", ".msh file name");
    parser.addPositionalArgument("solFileName", ".sol file name to render");
    parser.process(app);
    const QStringList args = parser.positionalArguments();
    if (args.size() != 2)
        {
        std::cerr << "error: expected a mesh file and a .sol file\n";
        parser.showHelp(1);
        }

    bool ok1, ok2;
    const double percent = parser.value("percent").toDouble(&ok1);
    const double qLength = parser.value("qlength").toDouble(&ok2);
    if (!ok1 || !(percent > 0.0) || percent > 100)
        {
        std::cerr << "error: -P percentage must be in ]0;100]\n";
        return 1;
        }
    if (!ok2)
        {
        std::cerr << "error: invalid quiver length\n";
        return 1;
        }
    const QString ext = "." + QFileInfo(args[0]).suffix().toLower();
    if (ext != ".msh")
        {
        std::cerr << ext.toStdString() << "  unsupported.\n";
        return 1;
        }

    std::array<int, 4> idx{};
    const bool subspace = parser.isSet("subspace");
    if (subspace)
        {
        const QStringList s =
                parser.value("subspace").split(' ', Qt::SkipEmptyParts);
        bool ok = s.size() == 4;
        for (int k = 0; ok && k < 4; ++k)
            {
            idx[k] = s[k].toInt(&ok);
            ok = ok && idx[k] >= 0 && idx[k] < 6;
            }
        if (!ok)
            {
            std::cerr << "error, --subspace option expects 4 indices in [0, 5]\n";
            return 1;
            }
        }

    std::vector<Node> nodes;
    try
        {
        fg::GmshSession gmshSession;
        const fg::Mesh mesh = fg::readMesh(args[0].toStdString(), true);
        const fg::TextTable sol = fg::readTextTable(args[1].toStdString());
        if (sol.rows.size() != mesh.points.size())
            {
            throw std::runtime_error("the .sol file has " + std::to_string(sol.rows.size())
                                     + " data lines, the mesh has "
                                     + std::to_string(mesh.points.size()) + " nodes");
            }
        for (std::size_t i = 0; i < mesh.points.size(); ++i)
            {
            const auto &p = mesh.points[i];
            const auto &r = sol.rows[i];
            if (r.size() < 4)
                {
                throw std::runtime_error("the .sol file has less than four columns");
                }
            nodes.push_back({p[0], p[1], p[2], r[1], r[2], r[3]});
            }
        nodes = filtering(nodes, parser.value("filter"), percent);
        }
    catch (const std::exception &e)
        {
        std::cerr << "error: " << e.what() << '\n';
        return 1;
        }

    QWidget *window = subspace ? plot2D(nodes, idx) : plot3D(nodes, qLength);
    window->setAttribute(Qt::WA_DeleteOnClose);
    window->setWindowTitle("fg-plot - " + QFileInfo(args[1]).fileName());
    window->show();
    return app.exec();
    }
