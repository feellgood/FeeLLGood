PROJECT (FEELLGOOD)
cmake_minimum_required(VERSION 2.8)

#===========================================================================
# FEELLGOOD version number.
#===========================================================================

set(FEELLGOOD_MAJOR_VERSION 0)
set(FEELLGOOD_MINOR_VERSION 3)
set(FEELLGOOD_PATCH_VERSION 0)

ENABLE_LANGUAGE(CXX)

MESSAGE(STATUS "CXX ${CMAKE_CXX_COMPILER_ID}")

#===========================================================================
# Path
#===========================================================================

SET(EXPRTK_INCLUDE_DIR ../exprtk )

SET(ANN_INCLUDE_DIR ../ann_1.1.2/include/ANN )

SET(GMM_INCLUDE_DIR ../gmm-5.3/include)

SET(SCALFMM_INCLUDE_DIR ../SCALFMM-1.4-148/Src )

SET(SCALFMM_CFG /usr/local/include/Utils )

INCLUDE_DIRECTORIES(${ANN_INCLUDE_DIR} ${GMM_INCLUDE_DIR} ${SCALFMM_INCLUDE_DIR} ${SCALFMM_CFG} ${EXPRTK_INCLUDE_DIR} )

#===========================================================================
# these parameters are used to build a config.h file 
#===========================================================================

set(FEELLGOOD_VERSION "\"${FEELLGOOD_MAJOR_VERSION}.${FEELLGOOD_MINOR_VERSION}.${FEELLGOOD_PATCH_VERSION}\"") 

set(my_cmd "git")
set(my_arg "show-ref")
set(my_arg1 "--head")
set(my_arg2 "HEAD")
set(my_arg3 "--hash")

execute_process(COMMAND ${my_cmd} ${my_arg} ${my_arg1} ${my_arg2} ${my_arg3} RESULT_VARIABLE git_result OUTPUT_VARIABLE git_var)

string(REPLACE "\n" "" SHAstr "\"${git_var}\"") 
MESSAGE(STATUS "version = ${FEELLGOOD_VERSION}")
MESSAGE(STATUS "executable built with SHA1=${SHAstr}")

set(VERBOSE_MODE true)

set(my_cmd "getconf")
set(my_arg "_NPROCESSORS_ONLN")
execute_process(COMMAND ${my_cmd} ${my_arg} RESULT_VARIABLE getnbcpu_result OUTPUT_VARIABLE MAX_NB_THREADS)
MESSAGE(STATUS "MAX_NB_THREADS=${MAX_NB_THREADS}")

#===========================================================================
# Compiler options
#===========================================================================

set( CMAKE_CXX_FLAGS_DEBUG "-DNDEBUG -Wall -Wextra -std=c++11 -fopenmp -O2 -ffast-math -msse4.2 -mfpmath=sse" CACHE STRING "C++ compiler flags used during debug builds." FORCE )

set( CMAKE_CXX_FLAGS_RELEASE "-DNDEBUG -Wall -Wextra -std=c++11 -fopenmp -O3 -ffast-math -msse4.2 -mfpmath=sse" CACHE STRING "C++ compiler flags used during release builds." FORCE )

set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -fopenmp")

SET( SOURCES config.h pt3D.h matBlocDiag.h tetra.h tetra.cpp facette.h facette.cpp feellgoodSettings.h node.h
	linear_algebra.h main.cc feellgoodSettings.cpp read.cpp femutil.cc
	recentrage.cc energy.cpp save.cpp vsolve.cc tiny.h)

SET( HEADERS ${ANN_INCLUDE_DIR}/ANN.h ${SCALFMM_INCLUDE_DIR}/Kernels/Rotation/FRotationKernel.hpp ${EXPRTK_INCLUDE_DIR}/exprtk.hpp )


configure_file(config.h.in ./config.h)

ADD_EXECUTABLE(feellgood ${SOURCES} ${HEADERS})

TARGET_LINK_LIBRARIES(feellgood scalfmm ANN)

install(TARGETS feellgood RUNTIME DESTINATION /usr/local/bin )

