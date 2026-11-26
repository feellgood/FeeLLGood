/** \file main.cpp
 * fg-3D: 3D viewer for gmsh mesh files and feeLLGood .sol files.
 */

#include <cstdlib>
#include <iostream>

#include <QApplication>
#include <QCommandLineParser>
#include <QSurfaceFormat>

#include <QVTKOpenGLNativeWidget.h>
#include "gmshMesh.h"

#include "mainWindow.h"

namespace
    {
/** True if nvidia-smi runs successfully, i.e. there is an NVIDIA GPU with its proprietary
 * driver. */
bool hasNvidiaGpu()
    {
#ifdef _WIN32
    return false;
#else
    return std::system("nvidia-smi > /dev/null 2>&1") == 0;
#endif
    }
    }  // namespace

/** Parse the command line and launch the GUI. On Unix-like systems, exit if there is no
 * display, and set the environment variables for NVIDIA PRIME render offload if an NVIDIA GPU
 * is detected. gmsh is initialized for the lifetime of the main window. */
int main(int argc, char *argv[])
    {
#ifndef _WIN32
    // Qt 5 uses X11 by default, unless QT_QPA_PLATFORM selects another platform
    const char *display = std::getenv("DISPLAY");
    const char *platform = std::getenv("QT_QPA_PLATFORM");
    if ((display == nullptr || *display == '\0') && (platform == nullptr || *platform == '\0'))
        {
        std::cerr << "Error: no display.\n";
        return 1;
        }
    if (hasNvidiaGpu())
        {
        setenv("__NV_PRIME_RENDER_OFFLOAD", "1", 1);
        setenv("__GLX_VENDOR_LIBRARY_NAME", "nvidia", 1);
        }
    else
        {
        std::cerr << "Warning: no NVIDIA GPU detected.\n";
        }
#endif

    QSurfaceFormat::setDefaultFormat(QVTKOpenGLNativeWidget::defaultFormat());
    QApplication app(argc, argv);
    QApplication::setApplicationName("fg-3D");
    QApplication::setApplicationVersion(FG_VERSION);

    QCommandLineParser parser;
    parser.setApplicationDescription(
            "fg-3D: 3D viewer for mesh and sol files\n\n"
            "Launch a graphic user interface to render magnetization from a feeLLGood .sol "
            "file as a vector field onto the corresponding mesh.");
    parser.addHelpOption();
    parser.addVersionOption();
    parser.addPositionalArgument(
            "meshFileName",
            "Optional gmsh mesh input file. If not provided, open one from the File menu.",
            "[meshFileName]");
    parser.process(app);
    const QStringList args = parser.positionalArguments();

    fg::GmshSession gmshSession;
    MainWindow window(args.isEmpty() ? QString() : args.first());
    window.show();
    return app.exec();
    }
