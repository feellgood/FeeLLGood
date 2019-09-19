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
# include directory path
#===========================================================================

SET(SCALFMM_CFG /usr/local/include/Utils )

INCLUDE_DIRECTORIES(${SCALFMM_CFG})

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

set(my_cmd "getconf")
set(my_arg "_NPROCESSORS_ONLN")
execute_process(COMMAND ${my_cmd} ${my_arg} RESULT_VARIABLE getnbcpu_result OUTPUT_VARIABLE MAX_NB_THREADS)
MESSAGE(STATUS "MAX_NB_THREADS=${MAX_NB_THREADS}")

#===========================================================================
# Compiler options
# -msse4.2 -mfpmath=sse (on X86 -mfpmath=sse is set by default ) or -mavx
#===========================================================================

set( CMAKE_CXX_FLAGS "-DNDEBUG -Wall -Wextra -std=c++11 -fopenmp -O3 -ffast-math -march=native" CACHE STRING "C++ compiler flags used during builds." FORCE )

set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -fopenmp")

SET( SOURCES config.h pt3D.h matBlocDiag.h tetra.h tetra.cpp facette.h facette.cpp feellgoodSettings.h node.h
	linear_algebra.h linear_algebra.cpp main.cpp feellgoodSettings.cpp read.cpp femutil.cpp
	recentering.cpp energy.cpp save.cpp solver.cpp tiny.h time_integration.cpp)


configure_file(config.h.in ./config.h)

ADD_EXECUTABLE(feellgood ${SOURCES})

TARGET_LINK_LIBRARIES(feellgood scalfmm ANN)

install(TARGETS feellgood RUNTIME DESTINATION /usr/local/bin )

