cmake_minimum_required(VERSION 3.13)
PROJECT(feellgood CXX)

# settings to add -std=c++17 option, and avoid -std=gnu++17
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

set(CMAKE_EXPORT_COMPILE_COMMANDS OFF) # ON to get a .json with compiler calls for all targets
# verbose might be done with make VERBOSE=1

find_package(OpenMP REQUIRED)

#===========================================================================
# cmake . -DENABLE_SEQ=OFF
# execution policy, by default it is parallel, sequential is for debug purpose.
# when ENABLE_SEQ is ON, executable does not link to TBB
#===========================================================================

option(ENABLE_SEQ "Enable sequential execution policy" OFF)
message( STATUS "Enable sequential execution policy: ENABLE_SEQ = " ${ENABLE_SEQ} )

if(NOT ${ENABLE_SEQ})
    find_package(TBB QUIET)
    if(TBB_FOUND)
        message(STATUS "Found TBB, parallel mode is active")
        set(TBB_LINK TBB::tbb)
    else()
        message(STATUS "TBB NOT FOUND")
        message(WARNING "sequential mode is active")
        set(ENABLE_SEQ TRUE)
    endif()
endif()

# Use the Intel MKL library if we can find it.
# First, try find_package(), which should work on a oneAPI environment.
set(MKL_INTERFACE lp64)
set(MKL_THREADING intel_thread)
find_package(MKL QUIET)
if(MKL_FOUND)
    message(STATUS "Found MKL with find_package(), MKL_CFLAGS= " ${MKL_CFLAGS})
    set(MKL_TARGET "MKL::MKL")
else()
    # As find_package() failed to find MKL, try pkg-config.
    find_package(PkgConfig)
    if(PKG_CONFIG_FOUND)
        pkg_check_modules(MKL QUIET IMPORTED_TARGET mkl-dynamic-lp64-iomp)
        if(MKL_FOUND)
            message(STATUS "Found MKL with pkg-config, MKL_CFLAGS= " ${MKL_CFLAGS})
            set(MKL_TARGET "PkgConfig::MKL")
        endif()
    else()
        message(STATUS "Not found: PkgConfig, needed for MKL")
    endif()
endif()

if(MKL_FOUND)
    set(EIGEN_MKL_CFLAGS "-DEIGEN_USE_MKL_ALL -DMKL_DIRECT_CALL")  # -DMKL_DIRECT_CALL : not sure we need it, neither it does something
else()
    message(STATUS "Not found: MKL")
endif()

find_package(Eigen3 3.3 REQUIRED NO_MODULE)

#===========================================================================
# these parameters are used to build a config.h file 
#===========================================================================

execute_process(
    COMMAND git describe --tags --always --dirty
    OUTPUT_VARIABLE FEELLGOOD_VERSION OUTPUT_STRIP_TRAILING_WHITESPACE
)
# remove leading "v" from tag name
string(REGEX REPLACE "^v" "" FEELLGOOD_VERSION ${FEELLGOOD_VERSION})

message(STATUS "feeLLGood version = ${FEELLGOOD_VERSION}")

# Link-time optimization requires cmake version at least 3.9.
option(ENABLE_LTO "Enable link-time optimization" OFF)
message(STATUS "Link-time optimization: ENABLE_LTO = " ${ENABLE_LTO})
if(${ENABLE_LTO})
    cmake_policy(SET CMP0069 NEW)
    include(CheckIPOSupported)
    check_ipo_supported()
    set(CMAKE_INTERPROCEDURAL_OPTIMIZATION TRUE)
endif()

set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -fopenmp")

SET(HEADERS config.h node.h expression_parser.h mesh.h electrostatSolver.h
    spinTransferTorque.h time_integration.h feellgoodSettings.h tetra.h
    facette.h triangle.h surface.h linear_algebra.h log-stats.h tags.h
    chronometer.h element.h)

SET(SOURCES feellgoodSettings.cpp time_integration.cpp solver.cpp
    read.cpp save.cpp linear_algebra.cpp recentering.cpp tetra.cpp
    energy.cpp facette.cpp expression_parser.cpp chronometer.cpp
    tags.cpp mesh.cpp)

configure_file(config.h.in ./config.h)

add_custom_command(
    OUTPUT default-settings.o
    COMMAND ld -r -b binary -z noexecstack default-settings.yml -o default-settings.o
    DEPENDS default-settings.yml)

#===========================================================================
# Compiler options
#===========================================================================
# -ffast-math -msse4.2 -mfpmath=sse (on X86 -mfpmath=sse set by default) or -mavx
# -fsanitize=address leak null undefined # both compiler and linker
#===========================================================================

ADD_EXECUTABLE(feellgood main.cpp ${HEADERS} ${SOURCES} default-settings.o)
target_compile_options( feellgood PUBLIC -Wall -Wextra -O3 -march=native ${MKL_CFLAGS} )
target_compile_definitions(feellgood PUBLIC -DNDEBUG -DTBB_SUPPRESS_DEPRECATED_MESSAGES ${EIGEN_MKL_CFLAGS} )

if (MKL_FOUND)
    #option --no-as-needed for mkl, see https://www.intel.com/content/www/us/en/developer/tools/oneapi/onemkl-link-line-advisor.html
    target_link_options(feellgood PUBLIC "-Wl,--no-as-needed" )
    # not sure we really need it since --no-as-needed is default value, though there might be a --as-needed before
    TARGET_LINK_LIBRARIES(feellgood yaml-cpp scalfmm duktape ANN OpenMP::OpenMP_CXX ${TBB_LINK} Eigen3::Eigen ${MKL_TARGET})
else()
    TARGET_LINK_LIBRARIES(feellgood yaml-cpp scalfmm duktape ANN OpenMP::OpenMP_CXX ${TBB_LINK} Eigen3::Eigen)
endif()

#===========================================================================
# Find where to install the Python modules.
#===========================================================================

# Find the version of the Python interpreter.
# Keep only the major and minor version numbers, e.g. "3.10".
execute_process(COMMAND python3 --version OUTPUT_VARIABLE PYTHON_VERSION)
string(
    REGEX REPLACE "Python ([0-9]+\\.[0-9]+)\\..*" "\\1"
    PYTHON_VERSION "${PYTHON_VERSION}"
)

# Possible install locations. The fisrt one (dist-packages) is used by the
# system Python on Debian and Debian-like systems. This is our preferred
# install location, and it will be used if it exists. The second one should
# work with any non-Debianized Python interpreter.
set(PY_DIST_DIR "lib/python${PYTHON_VERSION}/dist-packages")
set(PY_SITE_DIR "lib/python${PYTHON_VERSION}/site-packages")

# Select the best available installation directory.
if(EXISTS "${CMAKE_INSTALL_PREFIX}/${PY_DIST_DIR}")
    set(PY_MODULES_DIR "${PY_DIST_DIR}/feellgood")
elseif(PYTHON_VERSION)
    set(PY_MODULES_DIR "${PY_SITE_DIR}/feellgood")
endif()

# Install the python modules if we found a suitable location.
if(DEFINED PY_MODULES_DIR)
    install(DIRECTORY python-modules/ DESTINATION ${PY_MODULES_DIR})
    message(STATUS "Python modules will be installed in ${CMAKE_INSTALL_PREFIX}/${PY_MODULES_DIR}")
else()
    message(WARNING "Python modules will not be installed.")
endif()

#===========================================================================

install(TARGETS feellgood RUNTIME DESTINATION bin)
install(PROGRAMS ./tools/convert2vtk DESTINATION bin)
install(PROGRAMS ./tools/convert2hdf5 DESTINATION bin)
install(PROGRAMS ./tools/h5extract DESTINATION bin)
install(PROGRAMS ./tools/fg-merge-sol DESTINATION bin)

#===========================================================================
# unit tests, by default it is not built.
# you have to switch it ON while cmaking to build the unit tests and to be able to call 'make test' afterward:
#  cmake . -DENABLE_UTESTS=ON
# NB: ENABLE_UTESTS value is written in CMakeCache.txt
#===========================================================================

option(ENABLE_UTESTS "Enable tests" OFF)
message( STATUS "enable unit tests: ENABLE_UTESTS = "  ${ENABLE_UTESTS} )

#===========================================================================
#  cmake . -DENABLE_DET_UTESTS=ON
# this option forces random numbers to be deterministic 
#===========================================================================

option(ENABLE_DET_UTESTS "Enable deterministic unit tests" OFF)
message( STATUS "enable deterministic unit tests: ENABLE_DET_UTESTS = "  ${ENABLE_DET_UTESTS} )

if(${ENABLE_DET_UTESTS})
    set(DET_UTESTS true )
else()
    set(DET_UTESTS false )
endif()

if(${ENABLE_UTESTS})
    configure_file(./unit-tests/ut_config.h.in ./unit-tests/ut_config.h)
    include(CTest)
    add_subdirectory(unit-tests)
endif()

