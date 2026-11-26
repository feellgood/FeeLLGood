#include "mainWindow.h"

#include <algorithm>
#include <limits>

#include <QAction>
#include <QCheckBox>
#include <QCloseEvent>
#include <QDir>
#include <QDoubleSpinBox>
#include <QFileDialog>
#include <QFileInfo>
#include <QFont>
#include <QFrame>
#include <QGroupBox>
#include <QHBoxLayout>
#include <QLabel>
#include <QLocale>
#include <QMenuBar>
#include <QMessageBox>
#include <QSplitter>
#include <QStatusBar>
#include <QTextEdit>
#include <QVBoxLayout>

#include <QVTKOpenGLNativeWidget.h>
#include <vtkActor.h>
#include <vtkArrowSource.h>
#include <vtkAxesActor.h>
#include <vtkCamera.h>
#include <vtkCellArray.h>
#include <vtkCellData.h>
#include <vtkColorTransferFunction.h>
#include <vtkCubeAxesActor.h>
#include <vtkDoubleArray.h>
#include <vtkGenericOpenGLRenderWindow.h>
#include <vtkGlyph3D.h>
#include <vtkLookupTable.h>
#include <vtkOrientationMarkerWidget.h>
#include <vtkPointData.h>
#include <vtkPoints.h>
#include <vtkPolyData.h>
#include <vtkPolyDataMapper.h>
#include <vtkPolyDataNormals.h>
#include <vtkProperty.h>
#include <vtkRenderWindowInteractor.h>
#include <vtkRenderer.h>
#include <vtkScalarBarActor.h>
#include <vtkTextProperty.h>
#include <vtkUnstructuredGrid.h>
#include <vtkXMLUnstructuredGridWriter.h>

namespace
    {
const char *strHelp =
        "Workflow:\n"
        "Open mesh file from file menu. Only gmsh files can be read, any version.\n"
        "Open sol file from file menu. Only output FeeLLGood .sol files can be read.\n"
        "Export vtu : create a vtu (paraview) file from last opened mesh and sol files.\n";

const double RED[3] = {1, 0, 0};
const double GREEN[3] = {0, 0.5, 0};
const double BLUE[3] = {0, 0, 1};
const double BLACK[3] = {0, 0, 0};
const double LIGHTBLUE[3] = {0.678, 0.847, 0.902};
const double GREY[3] = {0.5, 0.5, 0.5};

/** Diverging blue-white-red lookup table, the 'coolwarm' colormap of matplotlib. */
vtkSmartPointer<vtkLookupTable> coolwarm(double vmin, double vmax)
    {
    vtkNew<vtkColorTransferFunction> ctf;
    ctf->SetColorSpaceToDiverging();
    ctf->AddRGBPoint(0.0, 0.230, 0.299, 0.754);
    ctf->AddRGBPoint(1.0, 0.706, 0.016, 0.150);
    const int n = 256;
    auto lut = vtkSmartPointer<vtkLookupTable>::New();
    lut->SetNumberOfTableValues(n);
    for (int i = 0; i < n; ++i)
        {
        double rgb[3];
        ctf->GetColor(static_cast<double>(i) / (n - 1), rgb);
        lut->SetTableValue(i, rgb[0], rgb[1], rgb[2], 1.0);
        }
    lut->SetTableRange(vmin, vmax);
    lut->Build();
    return lut;
    }
    }  // namespace

MainWindow::MainWindow(const QString &meshFile, QWidget *parent) : QMainWindow(parent)
    {
    setWindowTitle("FeeLLGood - 3D");
    setupUi();
    showMaximized();
    if (!meshFile.isEmpty())
        {
        openMesh(meshFile);
        }
    }

MainWindow::~MainWindow() = default;

void MainWindow::setupUi()
    {
    auto *central = new QWidget;
    setCentralWidget(central);
    strLog += "Welcome to fg-3D\n";
    auto *mainSplitter = new QSplitter(Qt::Horizontal);
    auto *centralLayout = new QHBoxLayout(central);
    centralLayout->addWidget(mainSplitter);
    centralLayout->setContentsMargins(2, 2, 2, 2);
    const int minLeftPanelWidth = 333;
    const int spacing = 2;

    // menus
    QMenu *fileMenu = menuBar()->addMenu("&File");
    QMenu *helpMenu = menuBar()->addMenu("&Help");

    QAction *meshAction = fileMenu->addAction("&Open .msh", this, &MainWindow::selectMshFile);
    meshAction->setStatusTip("Open a mesh file");
    meshAction->setShortcut(QKeySequence("Ctrl+M"));

    QAction *solAction = fileMenu->addAction("&Open .sol", this, &MainWindow::selectSolFile);
    solAction->setStatusTip("Open a simulation .sol file");
    solAction->setShortcut(QKeySequence("Ctrl+O"));

    QAction *helpAction = helpMenu->addAction("&Help", this, &MainWindow::showHelp);
    helpAction->setShortcut(QKeySequence("Ctrl+H"));

    fileMenu->addSeparator();

    QAction *exportAction = fileMenu->addAction("&Export vtu", this, &MainWindow::exportToVtu);
    exportAction->setStatusTip("Export 3D to vtu (paraview)");
    exportAction->setShortcut(QKeySequence("Ctrl+S"));

    // left panel
    auto *leftPanel = new QFrame;
    leftPanel->setFrameShape(QFrame::StyledPanel);
    leftPanel->setMinimumWidth(minLeftPanelWidth);
    auto *leftLayout = new QVBoxLayout(leftPanel);

    auto *rendererGroup = new QGroupBox("Renderer");
    auto *rendererLayout = new QVBoxLayout(rendererGroup);
    solLabel = new QLabel(QString::fromUtf8("✗ no .sol file selected."));
    rendererLayout->addWidget(solLabel);

    showVectorsCheckbox = new QCheckBox("Magnetization vectors");
    showVectorsCheckbox->setToolTip("add magnetization as arrows on the surface of the object");
    connect(showVectorsCheckbox, &QCheckBox::stateChanged, this, &MainWindow::updateView);

    surfaceChargesCheckbox = new QCheckBox("Surface charges");
    surfaceChargesCheckbox->setToolTip("add surface charges on the surface of the object");
    connect(surfaceChargesCheckbox, &QCheckBox::stateChanged, this,
            &MainWindow::onSurfaceChargesChanged);

    smoothShadingCheckbox = new QCheckBox("Smooth shading");
    smoothShadingCheckbox->setChecked(true);
    smoothShadingCheckbox->setToolTip("smooth the 3D rendering of the mesh");
    connect(smoothShadingCheckbox, &QCheckBox::stateChanged, this, &MainWindow::updateView);

    rendererLayout->addWidget(showVectorsCheckbox);
    rendererLayout->addWidget(surfaceChargesCheckbox);
    rendererLayout->addWidget(smoothShadingCheckbox);

    auto *scaleLayout = new QHBoxLayout;
    scaleLayout->addWidget(new QLabel("Vector Scale:"));
    vectorScaleSpinbox = new QDoubleSpinBox;
    vectorScaleSpinbox->setRange(0.001, 100.0);
    vectorScaleSpinbox->setValue(5.0);
    vectorScaleSpinbox->setSingleStep(0.1);
    vectorScaleSpinbox->setDecimals(3);
    vectorScaleSpinbox->setLocale(QLocale(QLocale::English, QLocale::UnitedStates));
    connect(vectorScaleSpinbox, &QDoubleSpinBox::editingFinished, this,
            &MainWindow::onVectorScaleChanged);
    scaleLayout->addWidget(vectorScaleSpinbox);
    rendererLayout->addLayout(scaleLayout);
    rendererLayout->setSpacing(spacing);
    leftLayout->addWidget(rendererGroup);

    auto *infosGroup = new QGroupBox("Renderer informations");
    auto *infosLayout = new QVBoxLayout(infosGroup);
    infosDisplay = new QTextEdit;
    infosDisplay->setToolTip("Informations on the current .msh and .sol files");
    infosDisplay->setReadOnly(true);
    infosLayout->addWidget(infosDisplay);
    infosLayout->setSpacing(spacing);
    leftLayout->addWidget(infosGroup, 1);
    leftLayout->setContentsMargins(2, 2, 2, 2);
    leftLayout->setSpacing(spacing);

    // right panel: 3D view
    vtkWidget = new QVTKOpenGLNativeWidget;
    vtkWidget->setRenderWindow(renderWindow);
    renderWindow->AddRenderer(renderer);
    renderer->SetBackground(1, 1, 1);
    renderer->GetActiveCamera()->ParallelProjectionOn();

    vtkNew<vtkAxesActor> axes;
    axesWidget->SetOrientationMarker(axes);
    axesWidget->SetInteractor(renderWindow->GetInteractor());
    axesWidget->SetViewport(0.0, 0.0, 0.2, 0.2);
    axesWidget->SetEnabled(1);
    axesWidget->InteractiveOff();

    mainSplitter->addWidget(leftPanel);
    mainSplitter->addWidget(vtkWidget);
    mainSplitter->setSizes({minLeftPanelWidth, 900});

    statusBar()->setFont(QFont("Times", 12));
    infosDisplay->setPlainText(strLog);
    updateView();
    }

void MainWindow::log(const QString &msg, bool status)
    {
    strLog += msg + "\n";
    infosDisplay->setPlainText(strLog);
    if (status)
        {
        statusBar()->showMessage(msg);
        }
    }

void MainWindow::closeEvent(QCloseEvent *event)
    {
    if (helpWindow)
        {
        helpWindow->close();
        }
    event->accept();
    }

void MainWindow::showHelp()
    {
    if (!helpWindow)
        {
        helpWindow = new QWidget;
        helpWindow->setAttribute(Qt::WA_DeleteOnClose);
        auto *layout = new QVBoxLayout(helpWindow);
        layout->addWidget(new QLabel(strHelp));
        helpWindow->setWindowTitle(QString("fg-3D v") + FG_VERSION + " help");
        }
    helpWindow->show();
    helpWindow->raise();
    }

bool MainWindow::openMesh(const QString &path)
    {
    QFileInfo info(path);
    QDir::setCurrent(info.absolutePath());
    meshFileName = info.absoluteFilePath();
    try
        {
        std::string msg = loadMesh(meshFileName.toStdString(), mesh);
        meshLoaded = true;
        magnetization.clear();
        solLoaded = false;
        surfaceCharges.clear();
        buildSurface();
        resetCamera = true;
        const QSignalBlocker b1(showVectorsCheckbox), b2(surfaceChargesCheckbox);
        showVectorsCheckbox->setChecked(false);
        surfaceChargesCheckbox->setChecked(false);
        log(QString::fromStdString(msg));
        updateView();
        return true;
        }
    catch (const std::exception &e)
        {
        log(QString::fromUtf8("✗ Error loading mesh: ") + QString::fromStdString(e.what()));
        return false;
        }
    }

void MainWindow::buildSurface()
    {
    vtkNew<vtkPoints> pts;
    pts->SetDataTypeToDouble();
    pts->SetNumberOfPoints(static_cast<vtkIdType>(mesh.points.size()));
    for (std::size_t i = 0; i < mesh.points.size(); ++i)
        {
        pts->SetPoint(static_cast<vtkIdType>(i), mesh.points[i].data());
        }
    vtkNew<vtkCellArray> tris;
    for (const auto &t : mesh.triangles)
        {
        vtkIdType ids[3] = {static_cast<vtkIdType>(t[0]), static_cast<vtkIdType>(t[1]),
                            static_cast<vtkIdType>(t[2])};
        tris->InsertNextCell(3, ids);
        }
    surface = vtkSmartPointer<vtkPolyData>::New();
    surface->SetPoints(pts);
    surface->SetPolys(tris);
    cellNormals = outwardNormals(mesh);
    }

void MainWindow::selectSolFile()
    {
    QString path = QFileDialog::getOpenFileName(this, "Select a .sol file", "",
                                                "SOL files (*.sol)");
    if (path.isEmpty())
        {
        return;
        }
    solFileName = path;
    solLabel->setText(QFileInfo(path).fileName());
    log(QString::fromUtf8("✓ .sol file ") + path + " selected");
    try
        {
        magnetization = loadSol(path.toStdString());
        solLoaded = true;
        surfaceCharges.clear();
        if (surfaceChargesCheckbox->isChecked())
            {
            onSurfaceChargesChanged();
            }
        else
            {
            updateView();
            }
        }
    catch (const std::exception &e)
        {
        QString msg = QString("Error loading .sol file: ") + e.what();
        QMessageBox::critical(this, "Error", msg);
        log(msg, false);
        magnetization.clear();
        solLoaded = false;
        }
    }

void MainWindow::selectMshFile()
    {
    QString path = QFileDialog::getOpenFileName(this, "Select a .msh file", "",
                                                "mesh files (*.msh)");
    if (!path.isEmpty())
        {
        openMesh(path);
        }
    }

void MainWindow::exportToVtu()
    {
    if (solFileName.isEmpty() || !solLoaded)
        {
        log(QString::fromUtf8(
                "✗ No magnetization loaded: open a .sol file before exporting to .vtu"));
        return;
        }
    if (!meshLoaded)
        {
        log(QString::fromUtf8("✗ No mesh loaded: open a .msh file before exporting to .vtu"));
        return;
        }
    const std::size_t nbNodes = mesh.points.size();
    vtkNew<vtkDoubleArray> mag;
    mag->SetName("magnetization");
    mag->SetNumberOfComponents(3);
    mag->SetNumberOfTuples(static_cast<vtkIdType>(nbNodes));
    for (std::size_t i = 0; i < nbNodes; ++i)
        {
        auto it = magnetization.find(i);
        if (it == magnetization.end())
            {
            log(QString::fromUtf8("✗ No magnetization for node %1 in the .sol file: "
                                  "the .sol file does not match the mesh")
                        .arg(i));
            return;
            }
        mag->SetTuple(static_cast<vtkIdType>(i), it->second.data());
        }

    vtkNew<vtkUnstructuredGrid> grid;
    grid->SetPoints(surface->GetPoints());
    grid->Allocate(static_cast<vtkIdType>(mesh.triangles.size() + mesh.tetrahedrons.size()));
    for (const auto &t : mesh.triangles)
        {
        vtkIdType ids[3] = {static_cast<vtkIdType>(t[0]), static_cast<vtkIdType>(t[1]),
                            static_cast<vtkIdType>(t[2])};
        grid->InsertNextCell(VTK_TRIANGLE, 3, ids);
        }
    for (const auto &t : mesh.tetrahedrons)
        {
        vtkIdType ids[4] = {static_cast<vtkIdType>(t[0]), static_cast<vtkIdType>(t[1]),
                            static_cast<vtkIdType>(t[2]), static_cast<vtkIdType>(t[3])};
        grid->InsertNextCell(VTK_TETRA, 4, ids);
        }
    grid->GetPointData()->AddArray(mag);

    QString outputFileName = QFileInfo(solFileName).completeBaseName() + ".vtu";
    vtkNew<vtkXMLUnstructuredGridWriter> writer;
    writer->SetFileName(outputFileName.toStdString().c_str());
    writer->SetInputData(grid);
    if (writer->Write() == 1)
        {
        log(QDir::current().filePath(outputFileName) + " written");
        }
    else
        {
        log(QString::fromUtf8("✗ Could not write ") + outputFileName);
        }
    }

void MainWindow::onSurfaceChargesChanged()
    {
    if (surfaceChargesCheckbox->isChecked())
        {
        if (meshLoaded && solLoaded)
            {
            surfaceCharges = computeSurfaceCharges(mesh, magnetization, cellNormals);
            if (surfaceCharges.empty())
                {
                log(QString::fromUtf8(
                        "✗ the .sol file does not match the mesh: no surface charges"));
                }
            else
                {
                log("surface charges computed", false);
                }
            }
        else
            {
            log("no mesh and/or magnetization data to compute surface charges", false);
            }
        }
    updateView();
    }

void MainWindow::onVectorScaleChanged()
    {
    if (showVectorsCheckbox->isChecked())
        {
        updateView();
        }
    }

void MainWindow::updateView()
    {
    renderer->RemoveAllViewProps();
    if (!meshLoaded || surface == nullptr)
        {
        statusBar()->showMessage(QString::fromUtf8("✗ Nothing to 3D render."));
        renderWindow->Render();
        return;
        }

    const bool smooth = smoothShadingCheckbox->isChecked();
    const bool showCharges = surfaceChargesCheckbox->isChecked() && !surfaceCharges.empty();

    vtkNew<vtkPolyData> displayed;
    displayed->ShallowCopy(surface);
    if (showCharges)
        {
        vtkNew<vtkDoubleArray> charges;
        charges->SetName("surface charges");
        charges->SetNumberOfTuples(static_cast<vtkIdType>(surfaceCharges.size()));
        for (std::size_t i = 0; i < surfaceCharges.size(); ++i)
            {
            charges->SetValue(static_cast<vtkIdType>(i), surfaceCharges[i]);
            }
        // a new cell data object, so that surface keeps its own
        vtkNew<vtkCellData> cd;
        cd->ShallowCopy(surface->GetCellData());
        cd->SetScalars(charges);
        displayed->GetCellData()->ShallowCopy(cd);
        }

    vtkNew<vtkPolyDataMapper> mapper;
    if (smooth)
        {  // point normals, sharp edges split
        vtkNew<vtkPolyDataNormals> pn;
        pn->SetInputData(displayed);
        pn->ComputePointNormalsOn();
        pn->ComputeCellNormalsOff();
        pn->SplittingOn();
        pn->SetFeatureAngle(30.0);
        pn->Update();
        mapper->SetInputConnection(pn->GetOutputPort());
        }
    else
        {
        mapper->SetInputData(displayed);
        }

    vtkNew<vtkActor> actor;
    actor->SetMapper(mapper);
    actor->GetProperty()->SetInterpolation(smooth ? VTK_PHONG : VTK_FLAT);
    if (showCharges)
        {
        auto [mn, mx] = std::minmax_element(surfaceCharges.begin(), surfaceCharges.end());
        auto lut = coolwarm(*mn, *mx);
        mapper->SetLookupTable(lut);
        mapper->SetScalarModeToUseCellData();
        mapper->SetScalarRange(*mn, *mx);
        mapper->ScalarVisibilityOn();

        vtkNew<vtkScalarBarActor> bar;
        bar->SetLookupTable(lut);
        bar->SetTitle("surface charges");
        bar->SetNumberOfLabels(5);
        bar->SetOrientationToHorizontal();
        bar->SetPosition(0.3, 0.02);
        bar->SetWidth(0.4);
        bar->SetHeight(0.08);
        bar->GetTitleTextProperty()->SetColor(BLACK[0], BLACK[1], BLACK[2]);
        bar->GetLabelTextProperty()->SetColor(BLACK[0], BLACK[1], BLACK[2]);
        renderer->AddActor2D(bar);
        }
    else
        {
        mapper->ScalarVisibilityOff();
        actor->GetProperty()->SetColor(LIGHTBLUE[0], LIGHTBLUE[1], LIGHTBLUE[2]);
        }
    renderer->AddActor(actor);

    if (!showCharges)
        {
        addBoundingBox();
        }
    if (showVectorsCheckbox->isChecked() && solLoaded)
        {
        addMagnetizationVectors();
        }
    addCubeAxes();

    if (resetCamera)
        {
        renderer->ResetCamera();
        resetCamera = false;
        }
    renderWindow->Render();
    }

void MainWindow::addMagnetizationVectors()
    {
    const double scale = vectorScaleSpinbox->value();
    vtkNew<vtkPoints> pts;
    pts->SetDataTypeToDouble();
    vtkNew<vtkDoubleArray> vectors;
    vectors->SetName("vectors");
    vectors->SetNumberOfComponents(3);
    for (std::size_t node : mesh.surfaceNodes)
        {
        auto it = magnetization.find(node);
        if (it == magnetization.end())
            {
            continue;
            }
        const Vec3 &p = mesh.points[node];
        const Vec3 &m = it->second;
        // arrows centered on the nodes
        pts->InsertNextPoint(p[0] - 0.5 * m[0] * scale, p[1] - 0.5 * m[1] * scale,
                             p[2] - 0.5 * m[2] * scale);
        double v[3] = {m[0] * scale, m[1] * scale, m[2] * scale};
        vectors->InsertNextTuple(v);
        }
    if (pts->GetNumberOfPoints() == 0)
        {
        statusBar()->showMessage(
                QString::fromUtf8("✗ No magnetization data found for surface nodes!"));
        return;
        }

    vtkNew<vtkPolyData> poly;
    poly->SetPoints(pts);
    poly->GetPointData()->SetVectors(vectors);

    vtkNew<vtkArrowSource> arrow;
    arrow->SetTipLength(0.25);
    arrow->SetTipRadius(0.1);
    arrow->SetShaftRadius(0.05);

    vtkNew<vtkGlyph3D> glyph;
    glyph->SetInputData(poly);
    glyph->SetSourceConnection(arrow->GetOutputPort());
    glyph->SetVectorModeToUseVector();
    glyph->OrientOn();
    glyph->SetScaleModeToDataScalingOff();
    glyph->SetScaleFactor(scale);

    vtkNew<vtkPolyDataMapper> mapper;
    mapper->SetInputConnection(glyph->GetOutputPort());
    mapper->ScalarVisibilityOff();
    vtkNew<vtkActor> actor;
    actor->SetMapper(mapper);
    actor->GetProperty()->SetColor(GREY[0], GREY[1], GREY[2]);
    actor->GetProperty()->SetOpacity(0.8);
    renderer->AddActor(actor);
    statusBar()->showMessage(QString::fromUtf8(
            "✓ Magnetization vectors plotted. If nothing is visible increase vector scale."));
    }

void MainWindow::addBoundingBox()
    {
    Vec3 lo, hi;
    lo.fill(std::numeric_limits<double>::max());
    hi.fill(std::numeric_limits<double>::lowest());
    for (const auto &p : mesh.points)
        {
        for (int k = 0; k < 3; ++k)
            {
            lo[k] = std::min(lo[k], p[k]);
            hi[k] = std::max(hi[k], p[k]);
            }
        }
    const Vec3 c[8] = {{lo[0], lo[1], lo[2]}, {hi[0], lo[1], lo[2]}, {hi[0], hi[1], lo[2]},
                       {lo[0], hi[1], lo[2]}, {lo[0], lo[1], hi[2]}, {hi[0], lo[1], hi[2]},
                       {hi[0], hi[1], hi[2]}, {lo[0], hi[1], hi[2]}};
    // each list holds three segments as pairs of corners
    addLines({c[0], c[1], c[3], c[2], c[4], c[5]}, RED, 4);
    addLines({c[1], c[2], c[0], c[3], c[4], c[7]}, GREEN, 4);
    addLines({c[0], c[4], c[1], c[5], c[3], c[7]}, BLUE, 4);
    addLines({c[7], c[6], c[5], c[6], c[2], c[6]}, BLACK, 2);
    }

void MainWindow::addLines(const std::vector<Vec3> &points, const double rgb[3], double width)
    {
    vtkNew<vtkPoints> pts;
    pts->SetDataTypeToDouble();
    vtkNew<vtkCellArray> lines;
    for (std::size_t i = 0; i + 1 < points.size(); i += 2)
        {
        vtkIdType ids[2] = {pts->InsertNextPoint(points[i].data()),
                            pts->InsertNextPoint(points[i + 1].data())};
        lines->InsertNextCell(2, ids);
        }
    vtkNew<vtkPolyData> poly;
    poly->SetPoints(pts);
    poly->SetLines(lines);
    vtkNew<vtkPolyDataMapper> mapper;
    mapper->SetInputData(poly);
    vtkNew<vtkActor> actor;
    actor->SetMapper(mapper);
    actor->GetProperty()->SetColor(rgb[0], rgb[1], rgb[2]);
    actor->GetProperty()->SetLineWidth(static_cast<float>(width));
    renderer->AddActor(actor);
    }

void MainWindow::addCubeAxes()
    {
    vtkNew<vtkCubeAxesActor> cube;
    cube->SetBounds(surface->GetBounds());
    cube->SetCamera(renderer->GetActiveCamera());
    cube->SetXTitle("X");
    cube->SetYTitle("Y");
    cube->SetZTitle("Z");
    cube->SetFlyModeToOuterEdges();
    cube->SetTickLocationToBoth();
    cube->SetGridLineLocation(vtkCubeAxesActor::VTK_GRID_LINES_FURTHEST);
    cube->DrawXGridlinesOn();
    cube->DrawYGridlinesOn();
    cube->DrawZGridlinesOn();
    cube->SetUseTextActor3D(0);
    for (int axis = 0; axis < 3; ++axis)
        {
        cube->GetTitleTextProperty(axis)->SetColor(BLACK[0], BLACK[1], BLACK[2]);
        cube->GetTitleTextProperty(axis)->SetFontSize(12);
        cube->GetLabelTextProperty(axis)->SetColor(BLACK[0], BLACK[1], BLACK[2]);
        cube->GetLabelTextProperty(axis)->SetFontSize(12);
        }
    const double grey[3] = {0.6, 0.6, 0.6};
    cube->GetXAxesLinesProperty()->SetColor(BLACK[0], BLACK[1], BLACK[2]);
    cube->GetYAxesLinesProperty()->SetColor(BLACK[0], BLACK[1], BLACK[2]);
    cube->GetZAxesLinesProperty()->SetColor(BLACK[0], BLACK[1], BLACK[2]);
    cube->GetXAxesGridlinesProperty()->SetColor(grey[0], grey[1], grey[2]);
    cube->GetYAxesGridlinesProperty()->SetColor(grey[0], grey[1], grey[2]);
    cube->GetZAxesGridlinesProperty()->SetColor(grey[0], grey[1], grey[2]);
    renderer->AddActor(cube);
    }
