#ifndef FG3D_MAINWINDOW_H
#define FG3D_MAINWINDOW_H

/** \file mainWindow.h
 * Main window of fg-3D: a control panel on the left, a VTK 3D view on the right.
 */

#include <QMainWindow>
#include <QPointer>
#include <QString>

#include <vtkNew.h>
#include <vtkSmartPointer.h>

#include "meshData.h"

class QCheckBox;
class QDoubleSpinBox;
class QLabel;
class QTextEdit;
class QVTKOpenGLNativeWidget;
class vtkActor;
class vtkGenericOpenGLRenderWindow;
class vtkOrientationMarkerWidget;
class vtkPolyData;
class vtkRenderer;

/** \class MainWindow
 * Main window of fg-3D. It holds the mesh and the magnetization currently loaded, and renders
 * them in a VTK view.
 */
class MainWindow : public QMainWindow
    {
    Q_OBJECT

public:
    /** Build the window, maximized. If meshFile is not empty, load and render it. */
    explicit MainWindow(const QString &meshFile = QString(), QWidget *parent = nullptr);

    /** Destructor, defined where the VTK classes are complete. */
    ~MainWindow() override;

protected:
    /** Close the help window, if any, along with the main window. */
    void closeEvent(QCloseEvent *event) override;

private slots:
    /** Slot of the "Open .msh" menu: ask for a mesh file, then load and render it. */
    void selectMshFile();

    /** Slot of the "Open .sol" menu: ask for a .sol file, load its magnetization and update
     * the view. Shows an error dialog if the file cannot be read. */
    void selectSolFile();

    /** Slot of the "Export vtu" menu: write the triangles, the tetrahedrons and the
     * magnetization on the nodes to a .vtu file named after the .sol file, in the directory
     * of the mesh file. Nothing is written if no .sol file is loaded or if it does not match
     * the mesh. */
    void exportToVtu();

    /** Slot of the "Help" menu: show the help window, creating it if needed. */
    void showHelp();

    /** Slot of the "Surface charges" checkbox: compute the surface charges when it is
     * checked and a mesh and a .sol file are loaded, then update the view. */
    void onSurfaceChargesChanged();

    /** Slot of the vector scale spinbox: update the view if the arrows are shown. */
    void onVectorScaleChanged();

    /** Render again the whole 3D view from the current data and options: surface, surface
     * charges or bounding box, magnetization arrows, grid and axes. The camera is reset only
     * after a new mesh is loaded. */
    void updateView();

private:
    /** Build the menus, the left control panel and the VTK view. */
    void setupUi();

    /** Load the mesh file at path, make its directory the current one, forget the previous
     * magnetization and render the mesh. Returns false and logs an error if it cannot be
     * loaded. */
    bool openMesh(const QString &path);

    /** Build the surface polydata from the triangles of the mesh, and store their outward
     * normals in cellNormals. */
    void buildSurface();

    /** Add the bounding box of the mesh: X edges in red, Y edges in green, Z edges in blue
     * from the minimum corner, the three edges of the maximum corner in black. */
    void addBoundingBox();

    /** Add line segments to the view. Consecutive pairs of points are the ends of a segment;
     * rgb is the color and width the line width in pixels. */
    void addLines(const std::vector<Vec3> &points, const double rgb[3], double width);

    /** Add the magnetization of the surface nodes as grey arrows, centered on the nodes and
     * scaled by the value of the vector scale spinbox. */
    void addMagnetizationVectors();

    /** Add the X, Y, Z axes with ticks, labels and grid lines on the far side of the mesh
     * bounding box. */
    void addCubeAxes();

    /** Append msg to the information text box and, if status is true, show it in the status
     * bar. */
    void log(const QString &msg, bool status = true);

    QString meshFileName;
    QString solFileName;
    MeshData mesh;
    bool meshLoaded = false;
    Magnetization magnetization;
    bool solLoaded = false;
    std::vector<Vec3> cellNormals;  //!< outward normals of the triangles
    std::vector<double> surfaceCharges;
    bool resetCamera = true;
    QString strLog;

    QLabel *solLabel = nullptr;
    QCheckBox *showVectorsCheckbox = nullptr;
    QCheckBox *surfaceChargesCheckbox = nullptr;
    QCheckBox *smoothShadingCheckbox = nullptr;
    QDoubleSpinBox *vectorScaleSpinbox = nullptr;
    QTextEdit *infosDisplay = nullptr;
    QPointer<QWidget> helpWindow;

    QVTKOpenGLNativeWidget *vtkWidget = nullptr;
    vtkNew<vtkGenericOpenGLRenderWindow> renderWindow;
    vtkNew<vtkRenderer> renderer;
    vtkNew<vtkOrientationMarkerWidget> axesWidget;
    vtkSmartPointer<vtkPolyData> surface;  //!< triangles of the mesh
    };

#endif
