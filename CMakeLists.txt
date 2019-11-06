PROJECT (FEELLGOOD)
cmake_minimum_required(VERSION 2.8)

#===========================================================================
# FEELLGOOD version number
#===========================================================================

set(FEELLGOOD_MAJOR_VERSION 0)
set(FEELLGOOD_MINOR_VERSION 6)
set(FEELLGOOD_PATCH_VERSION 0)

ENABLE_LANGUAGE(CXX)

MESSAGE(STATUS "CXX ${CMAKE_CXX_COMPILER_ID}")

#===========================================================================
# these parameters are used to build a config.h file 
#===========================================================================

set(FEELLGOOD_VERSION "\"${FEELLGOOD_MAJOR_VERSION}.${FEELLGOOD_MINOR_VERSION}.${FEELLGOOD_PATCH_VERSION}\"") 

execute_process(
    COMMAND git show-ref --verify --head HEAD --hash
    OUTPUT_VARIABLE SHAstr OUTPUT_STRIP_TRAILING_WHITESPACE
)
set(SHAstr "\"${SHAstr}\"")  # add quotes

MESSAGE(STATUS "version = ${FEELLGOOD_VERSION} from SHA1 = ${SHAstr}")

execute_process(
    COMMAND getconf _NPROCESSORS_ONLN
    OUTPUT_VARIABLE MAX_NB_THREADS OUTPUT_STRIP_TRAILING_WHITESPACE
)
MESSAGE(STATUS "MAX_NB_THREADS = ${MAX_NB_THREADS}")

#===========================================================================
# Compiler options -ffast-math
# -msse4.2 -mfpmath=sse (on X86 -mfpmath=sse is set by default ) or -mavx
#===========================================================================

set( CMAKE_CXX_FLAGS "-DNDEBUG -Wall -Wextra -std=c++11 -fopenmp -O3 -march=native -fsanitize=address" CACHE STRING "C++ compiler flags used during builds." FORCE )

set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -fopenmp")

SET(HEADERS config.h pt3D.h node.h matBlocDiag.h mag_parser.h
    time_integration.h feellgoodSettings.h tetra.h facette.h linear_algebra.h)

SET(SOURCES mag_parser.cpp time_integration.cpp feellgoodSettings.cpp
    solver.cpp tetra.cpp facette.cpp read.cpp save.cpp
    linear_algebra.cpp femutil.cpp energy.cpp recentering.cpp tiny.h)

configure_file(config.h.in ./config.h)

ADD_EXECUTABLE(feellgood main.cpp ${HEADERS} ${SOURCES})

TARGET_LINK_LIBRARIES(feellgood scalfmm ANN)

if(DEFINED withSharedLib)
    MESSAGE( "Shared lib generation " ${withSharedLib} )
    if (withSharedLib)
        ADD_LIBRARY( _feellgood SHARED ${HEADERS} ${SOURCES} )
        install(FILES feellgood_core.h DESTINATION include)
        install(TARGETS _feellgood DESTINATION lib)
    endif(withSharedLib)
endif(DEFINED withSharedLib)

install(TARGETS feellgood RUNTIME DESTINATION bin)

