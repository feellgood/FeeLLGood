SET(HEADERS config.h node.h expression_parser.h mesh.h electrostatSolver.h
        spinAccumulationSolver.h time_integration.h settings.h tetra.h
        facette.h linear_algebra.h log-stats.h tags.h chronometer.h element.h)

SET(SOURCES settings.cpp time_integration.cpp solver.cpp electrostatSolver.cpp
        spinAccumulationSolver.cpp read.cpp save.cpp linear_algebra.cpp recentering.cpp tetra.cpp
        energy.cpp facette.cpp expression_parser.cpp chronometer.cpp
        tags.cpp mesh.cpp)

configure_file(config.h.in ./config.h)

set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -fopenmp")

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

include_directories (${CMAKE_SOURCE_DIR} ${GMSH_INC})

ADD_EXECUTABLE(feellgood main.cpp ${HEADERS} ${SOURCES} default-settings.o)
target_compile_options( feellgood PUBLIC -Wall -Wextra -O3 -march=native ${MKL_CFLAGS} )
target_compile_definitions(feellgood PUBLIC
        -DNDEBUG -DTBB_SUPPRESS_DEPRECATED_MESSAGES ${EIGEN_MKL_CFLAGS} )

if (MKL_FOUND)
    #option --no-as-needed for mkl, see
    # https://www.intel.com/content/www/us/en/developer/tools/oneapi/onemkl-link-line-advisor.html
    target_link_options(feellgood PUBLIC "-Wl,--no-as-needed" )
    # not sure we really need it since --no-as-needed is default value, though there might be
    # a --as-needed before
    TARGET_LINK_LIBRARIES(feellgood yaml-cpp duktape ann OpenMP::OpenMP_CXX ${TBB_LINK}
            Eigen3::Eigen ${MKL_TARGET} ${GMSH_LIB})
else()
    TARGET_LINK_LIBRARIES(feellgood yaml-cpp duktape ann OpenMP::OpenMP_CXX ${TBB_LINK}
            Eigen3::Eigen ${GMSH_LIB})
endif()
