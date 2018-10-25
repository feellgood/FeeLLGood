PROJECT (FEELLGOOD)
cmake_minimum_required(VERSION 2.8)

#===========================================================================
# FEELLGOOD version number.
#===========================================================================

set(FEELLGOOD_MAJOR_VERSION 0)
set(FEELLGOOD_MINOR_VERSION 0)
set(FEELLGOOD_PATCH_VERSION 0)
set(FEELLGOOD_VERSION  "${FEELLGOOD_MAJOR_VERSION}.${FEELLGOOD_MINOR_VERSION}.${FEELLGOOD_PATCH_VERSION}" )

ENABLE_LANGUAGE(CXX)

MESSAGE(STATUS "CXX ${CMAKE_CXX_COMPILER_ID}")

#===========================================================================
# Path
#===========================================================================

SET(ANN_INCLUDE_DIR ../ann_1.1.2/include/ANN )

SET(MTL_INCLUDE_DIR ../mtl4 )

SET(SCALFMM_INCLUDE_DIR ../SCALFMM-1.4-148/Src )

SET(SCALFMM_CFG /usr/local/include/Utils )

INCLUDE_DIRECTORIES(${ANN_INCLUDE_DIR} ${GMM_INCLUDE_DIR} ${SCALFMM_INCLUDE_DIR} ${SCALFMM_CFG})

#===========================================================================
# these parameters are used to build a config.h file 
#===========================================================================

set(git_cmd "git")
set(git_arg "show-ref")
set(git_arg1 "--head")
set(git_arg2 "HEAD")
set(git_arg3 "--hash")

execute_process(COMMAND ${git_cmd} ${git_arg} ${git_arg1} ${git_arg2} ${git_arg3} RESULT_VARIABLE git_result OUTPUT_VARIABLE git_var)

string(REPLACE "\n" "" SHAstr "\"${git_var}\"") 
MESSAGE(STATUS "executable built with SHA1=${SHAstr}")


set(VERBOSE_MODE true)
set(SECOND_ORDER true)
set(NB_THREADS 16)
MESSAGE(STATUS "NB_THREADS=${NB_THREADS}")
#===========================================================================
# Compiler options
#===========================================================================

set( CMAKE_CXX_FLAGS_DEBUG "-Wall -Wextra -std=c++11 -O2 -fopenmp -ffast-math -msse4.2 -mfpmath=sse" CACHE STRING "C++ compiler flags used during debug builds." FORCE )

set( CMAKE_CXX_FLAGS_RELEASE "-Wall -Wextra -std=c++11 -O3 -fopenmp -ffast-math -msse4.2 -mfpmath=sse" CACHE STRING "C++ compiler flags used during release builds." FORCE )

SET( SOURCES config.h pt3D.h matBlocDiag.h tetra.h tetra.cpp facette.h facette.cpp feellgoodSettings.h node.h
	linear_algebra.h main.cc feellgoodSettings.cpp read.cpp femutil.cc
	recentrage.cc energy.cc save.cpp vsolve.cc tiny.h)

SET( HEADERS ${ANN_INCLUDE_DIR}/ANN.h ${SCALFMM_INCLUDE_DIR}/Kernels/Rotation/FRotationKernel.hpp)


configure_file(config.h.in ./config.h)

ADD_EXECUTABLE(feellgood ${SOURCES} ${HEADERS})

TARGET_LINK_LIBRARIES(feellgood scalfmm ANN)

install(TARGETS feellgood RUNTIME DESTINATION /usr/local/bin )

