cmake_minimum_required(VERSION 3.16)
project(feellgood-tools VERSION 1.0.0 LANGUAGES C CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_AUTOMOC ON)

if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release)
endif()

# The Qt version must match the one VTK's GUISupportQt module was built with.
find_package(Qt5 REQUIRED COMPONENTS Core Gui Widgets Qml)

find_library(GMSH_LIB gmsh)
if(NOT GMSH_LIB)
    message(FATAL_ERROR "Could not find libgmsh")
endif()
find_path(GMSH_INC gmsh.h)
if(NOT GMSH_INC)
    message(FATAL_ERROR "Could not find gmsh.h")
endif()

# Optional dependencies: the tools needing a missing one are not built.

# Some VTK packages (Debian, Ubuntu) link VTK::mpi to MPI::MPI_C without
# finding MPI, hence the C language in project().
find_package(MPI QUIET COMPONENTS C)
find_package(VTK QUIET COMPONENTS
    CommonCore
    CommonDataModel
    FiltersCore
    FiltersSources
    GUISupportQt
    InteractionStyle
    InteractionWidgets
    IOXML
    RenderingAnnotation
    RenderingCore
    RenderingFreeType
    RenderingOpenGL2)
if(NOT VTK_FOUND)
    message(WARNING "VTK not found: the tools needing it will not be built")
endif()

add_subdirectory(common)

# Add an executable named after its directory, with the version given.
# Extra arguments are libraries to link with, besides the common library.
function(add_fg_tool name version)
    file(GLOB sources ${CMAKE_CURRENT_SOURCE_DIR}/${name}/*.cpp)
    add_executable(${name} ${sources})
    target_compile_definitions(${name} PRIVATE FG_VERSION="${version}")
    target_link_libraries(${name} PRIVATE fgcommon ${ARGN})
    install(TARGETS ${name} RUNTIME DESTINATION bin)
endfunction()

add_fg_tool(fg-gensol 1.0.1)
add_fg_tool(fg-merge-sol 1.0.0)

if(VTK_FOUND)
    add_fg_tool(fg-3D 0.0.0 Qt5::Widgets ${VTK_LIBRARIES})
    vtk_module_autoinit(TARGETS fg-3D MODULES ${VTK_LIBRARIES})
    add_fg_tool(fg-plot 1.0.0 Qt5::Widgets ${VTK_LIBRARIES})
    vtk_module_autoinit(TARGETS fg-plot MODULES ${VTK_LIBRARIES})
    add_fg_tool(fg-sol2vtk 1.1.2 ${VTK_LIBRARIES})
    vtk_module_autoinit(TARGETS fg-sol2vtk MODULES ${VTK_LIBRARIES})
endif()
