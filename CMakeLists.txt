PROJECT (FEELLGOOD)
cmake_minimum_required(VERSION 2.8)

ENABLE_LANGUAGE(CXX)

MESSAGE(STATUS "CXX ${CMAKE_CXX_COMPILER_ID}")

SET(ANN_INCLUDE_DIR /home/nanospin/Documents/ann_1.1.2/include/ANN)

SET(GMM_INCLUDE_DIR /home/nanospin/Documents/gmm-5.2/include/gmm)

SET(SCALFMM_INCLUDE_DIR /home/nanospin/Documents/ScalFMM_1.5/Src)

INCLUDE_DIRECTORIES(${ANN_INCLUDE_DIR} ${GMM_INCLUDE_DIR} ${SCALFMM_INCLUDE_DIR})

SET(INCLUDE ${ANN_INCLUDE_DIR} ${GMM_INCLUDE_DIR} ${SCALFMM_INCLUDE_DIR})

set( CMAKE_CXX_FLAGS "-DORD2 -DNDEBUG -DNSTAT -Wall -frounding-math -std=c++11 -O3 -fopenmp -ffast-math -msse4.2 -mfpmath=sse" CACHE STRING "C++ compiler flags used during builds." FORCE )

set( CMAKE_CXX_FLAGS_DEBUG "-DORD1 -DNDEBUG -DNSTAT -Wall -frounding-math -std=c++11 -O2 -fopenmp -ffast-math -msse4.2 -mfpmath=sse" CACHE STRING "C++ compiler flags used during builds." FORCE )

set( CMAKE_CXX_FLAGS_RELEASE "-DORD2 -DNDEBUG -DNSTAT -Wall -frounding-math -std=c++11 -O3 -fopenmp -ffast-math -msse4.2 -mfpmath=sse" CACHE STRING "C++ compiler flags used during builds." FORCE )

# file(GLOB sources . *.cc)

SET( SOURCES main.cc init_distrib.cc cputime.cc extract_comment.cc dialog.cc affichage.cc
	saveH.cc lecture.cc femutil.cc chapeaux.cc restoresol.cc direction.cc
	recentrage.cc u_moy.cc v_moy.cc energy.cc saver.cc savesol.cc savecfg_vtk.cc
	vsolve.cc base_projection.cc integrales.cc projection.cc assemblage.cc
	evolution.cc reset.cc tiny.h)

SET( HEADERS ${ANN_INCLUDE_DIR}/ANN.h ${GMM_INCLUDE_DIR}/gmm.h ${SCALFMM_INCLUDE_DIR}/Kernels/Rotation/FRotationKernel.hpp)


ADD_EXECUTABLE(feellgood_exe ${SOURCES} ${HEADERS})
SET( SCALFMM "libscalfmm.a")

#SET(ANN_LINK_DIR /home/nanospin/Documents/ann_1.1.2/lib)
SET( ANN "/home/nanospin/Documents/ann_1.1.2/lib/libANN.a")

#LINK_DIRECTORIES(${ANN_LINK_DIR}) init_distrib = quoi ça ???

TARGET_LINK_LIBRARIES(feellgood_exe ${SCALFMM} ${ANN})

