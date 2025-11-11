PROJECT (feellgood CXX)
message(STATUS " *** unit test cmaking *** ")
find_package (Boost COMPONENTS system filesystem unit_test_framework REQUIRED)
find_package(OpenMP REQUIRED)
find_package(TBB REQUIRED)
find_package(Eigen3 3.3 REQUIRED NO_MODULE)

find_library(GMSH_LIB gmsh)
if(NOT GMSH_LIB)
    message(STATUS "Could not find libgmsh")
    set(GMSH_TARGET "GMSH::GMSH")
endif()

find_path(GMSH_INC gmsh.h)
if(NOT GMSH_INC)
  message(FATAL_ERROR "Could not find gmsh.h")
endif()

include_directories (${CMAKE_SOURCE_DIR} ${Boost_INCLUDE_DIRS} ${GMSH_INC})

set( BOOST_TEST_LOG_LEVEL "all" )

add_definitions ( -DBOOST_TEST_DYN_LINK -DTBB_SUPPRESS_DEPRECATED_MESSAGES )
#add_compile_options( -fsanitize=leak )
#add_link_options( -fsanitize=leak )

add_executable (test_ut_solver ut_solver.cpp)

add_executable (test_ut_algebra_bicg ut_algebra_bicg.cpp)

add_executable (test_ut_algebra ut_algebra.cpp)

add_executable (test_ut_OMP_solver ut_OMP_solver.cpp)

add_executable(test_ut_dukTapeParser ../expression_parser.cpp ut_dukTapeParser.cpp)

add_executable (test_ut_node ut_node.cpp)

add_executable (test_ut_index ut_index.cpp)

SET(SOURCES ../facette.cpp ../tetra.cpp ut_element.cpp)
add_executable (test_ut_element ${SOURCES})

SET(SOURCES ../facette.cpp ut_facette.cpp)
add_executable (test_ut_facette ${SOURCES})

SET(SOURCES ../tetra.cpp ut_tetra.cpp)
add_executable (test_ut_tetra ${SOURCES})

SET(SOURCES ../tetra.cpp ut_tet_lumping.cpp)
add_executable (test_ut_tet_lumping ${SOURCES})

SET(SOURCES ../tetra.cpp ut_tet_charges.cpp)
add_executable (test_ut_tet_charges ${SOURCES})

SET(SOURCES ../tetra.cpp ut_energy.cpp ../facette.cpp)
add_executable (test_ut_energy ${SOURCES})

SET(SOURCES ../tetra.cpp ut_anisotropy.cpp)
add_executable (test_ut_anisotropy ${SOURCES})

add_executable (test_ut_tiny ut_tiny.cpp)

set(SOURCES ../tetra.cpp ut_time_int.cpp)
add_executable (test_ut_time_int ${SOURCES})

add_executable (test_ut_log-stats ut_log-stats.cpp)

add_executable(test_ut_readMesh ut_readMesh.cpp)

target_link_libraries(test_ut_solver
    ${Boost_FILESYSTEM_LIBRARY}
    ${Boost_SYSTEM_LIBRARY}
    ${Boost_UNIT_TEST_FRAMEWORK_LIBRARY}
    Eigen3::Eigen
    )

target_link_libraries(test_ut_element
    ${Boost_FILESYSTEM_LIBRARY}
    ${Boost_SYSTEM_LIBRARY}
    ${Boost_UNIT_TEST_FRAMEWORK_LIBRARY}
    TBB::tbb
    )

target_link_libraries(test_ut_tetra
    ${Boost_FILESYSTEM_LIBRARY}
    ${Boost_SYSTEM_LIBRARY}
    ${Boost_UNIT_TEST_FRAMEWORK_LIBRARY}
    TBB::tbb
    )

target_link_libraries(test_ut_anisotropy
    ${Boost_FILESYSTEM_LIBRARY}
    ${Boost_SYSTEM_LIBRARY}
    ${Boost_UNIT_TEST_FRAMEWORK_LIBRARY}
    TBB::tbb
    )

target_link_libraries(test_ut_tet_lumping
    ${Boost_FILESYSTEM_LIBRARY}
    ${Boost_SYSTEM_LIBRARY}
    ${Boost_UNIT_TEST_FRAMEWORK_LIBRARY}
    TBB::tbb
    )

target_link_libraries(test_ut_tet_charges
    ${Boost_FILESYSTEM_LIBRARY}
    ${Boost_SYSTEM_LIBRARY}
    ${Boost_UNIT_TEST_FRAMEWORK_LIBRARY}
    TBB::tbb
    )

target_link_libraries(test_ut_energy
    ${Boost_FILESYSTEM_LIBRARY}
    ${Boost_SYSTEM_LIBRARY}
    ${Boost_UNIT_TEST_FRAMEWORK_LIBRARY}
    TBB::tbb
    )

if (MKL_FOUND)
    #target_compile_options(test_ut_solver PUBLIC -fsanitize=leak )
    #target_link_options(test_ut_solver PUBLIC -fsanitize=leak )
    target_link_libraries(test_ut_solver PkgConfig::MKL)

    #target_compile_options(test_ut_element PUBLIC -fsanitize=leak )
    #target_link_options(test_ut_element PUBLIC -fsanitize=leak ) # -> SEGFAULT
    target_link_libraries(test_ut_element PkgConfig::MKL)

    #target_compile_options(test_ut_tetra PUBLIC -fsanitize=leak )
    #target_link_options(test_ut_tetra PUBLIC -fsanitize=leak )
    target_link_libraries(test_ut_tetra PkgConfig::MKL)

    target_link_libraries(test_ut_tet_lumping PkgConfig::MKL)

    target_link_libraries(test_ut_tet_charges PkgConfig::MKL)

    target_link_libraries(test_ut_energy PkgConfig::MKL)

    target_link_libraries(test_ut_anisotropy PkgConfig::MKL)
endif()

target_link_libraries(test_ut_algebra_bicg
  ${Boost_FILESYSTEM_LIBRARY}
  ${Boost_SYSTEM_LIBRARY}
  ${Boost_UNIT_TEST_FRAMEWORK_LIBRARY}
  TBB::tbb
  )

target_link_libraries(test_ut_algebra
        ${Boost_FILESYSTEM_LIBRARY}
        ${Boost_SYSTEM_LIBRARY}
        ${Boost_UNIT_TEST_FRAMEWORK_LIBRARY}
        TBB::tbb
        )

#target_compile_options(test_ut_OMP_solver PUBLIC -fsanitize=leak )
#target_link_options(test_ut_OMP_solver PUBLIC -fsanitize=leak )
target_link_libraries(test_ut_OMP_solver
  ${Boost_FILESYSTEM_LIBRARY}
  ${Boost_SYSTEM_LIBRARY}
  ${Boost_UNIT_TEST_FRAMEWORK_LIBRARY}
  Eigen3::Eigen
  OpenMP::OpenMP_CXX
  )

target_link_libraries(test_ut_dukTapeParser
  ${Boost_FILESYSTEM_LIBRARY}
  ${Boost_SYSTEM_LIBRARY}
  ${Boost_UNIT_TEST_FRAMEWORK_LIBRARY}
  duktape
  Eigen3::Eigen
  )

#target_compile_options(test_ut_node PUBLIC -fsanitize=leak )
#target_link_options(test_ut_node PUBLIC -fsanitize=leak ) # -> random SEGFAULT
target_link_libraries(test_ut_node
  ${Boost_FILESYSTEM_LIBRARY}
  ${Boost_SYSTEM_LIBRARY}
  ${Boost_UNIT_TEST_FRAMEWORK_LIBRARY}
  )

target_compile_options(test_ut_index PUBLIC -Wno-deprecated-declarations )
# warning about deprecated code is silented here because it is the goal of one test in test_ut_index
# to check the replacement of a transform(...,bind2nd(...)) by a loop
# the deprecated code is found in the LLG solver of feellgood research version for LLG + spin accumulation (June 2025)
target_link_libraries(test_ut_index
  ${Boost_FILESYSTEM_LIBRARY}
  ${Boost_SYSTEM_LIBRARY}
  ${Boost_UNIT_TEST_FRAMEWORK_LIBRARY}
  )

#target_compile_options(test_ut_facette PUBLIC -fsanitize=leak )
#target_link_options(test_ut_facette PUBLIC -fsanitize=leak )
target_link_libraries(test_ut_facette
  ${Boost_FILESYSTEM_LIBRARY}
  ${Boost_SYSTEM_LIBRARY}
  ${Boost_UNIT_TEST_FRAMEWORK_LIBRARY}
  TBB::tbb
  )

#target_compile_options(test_ut_tiny PUBLIC -fsanitize=leak )
#target_link_options(test_ut_tiny PUBLIC -fsanitize=leak )
target_link_libraries(test_ut_tiny
  ${Boost_FILESYSTEM_LIBRARY}
  ${Boost_SYSTEM_LIBRARY}
  ${Boost_UNIT_TEST_FRAMEWORK_LIBRARY}
  )

#target_compile_options(test_ut_time_int PUBLIC -fsanitize=leak )
#target_link_options(test_ut_time_int PUBLIC -fsanitize=leak )
target_link_libraries(test_ut_time_int
  ${Boost_FILESYSTEM_LIBRARY}
  ${Boost_SYSTEM_LIBRARY}
  ${Boost_UNIT_TEST_FRAMEWORK_LIBRARY}
  TBB::tbb
  )

#target_compile_options(test_ut_log-stats PUBLIC -fsanitize=leak )
#target_link_options(test_ut_log-stats PUBLIC -fsanitize=leak )
target_link_libraries(test_ut_log-stats
  ${Boost_UNIT_TEST_FRAMEWORK_LIBRARY}
  )

target_link_libraries(test_ut_readMesh
  ${Boost_UNIT_TEST_FRAMEWORK_LIBRARY}
  ${GMSH_LIB}
  TBB::tbb
  )

# setarch --addr-no-randomize ./test_ut_* to disable ASLR to avoid -fsanitize=leak to bug
add_test (NAME ut_solver COMMAND test_ut_solver)
add_test (NAME ut_algebra_bicg COMMAND test_ut_algebra_bicg)
add_test (NAME ut_OMP_solver COMMAND test_ut_OMP_solver)
add_test (NAME ut_algebra COMMAND test_ut_algebra)
add_test (NAME ut_dukTapeParser COMMAND test_ut_dukTapeParser)
add_test (NAME ut_node COMMAND test_ut_node)
add_test (NAME ut_index COMMAND test_ut_index)
add_test (NAME ut_element COMMAND test_ut_element)
add_test (NAME ut_facette COMMAND test_ut_facette)
add_test (NAME ut_tetra COMMAND test_ut_tetra)
add_test (NAME ut_tet_lumping COMMAND test_ut_tet_lumping)
add_test (NAME ut_tet_charges COMMAND test_ut_tet_charges)
add_test (NAME ut_energy COMMAND test_ut_energy)
add_test (NAME ut_anisotropy COMMAND test_ut_anisotropy)
add_test (NAME ut_tiny COMMAND test_ut_tiny)
add_test (NAME ut_time_int COMMAND test_ut_time_int)
add_test (NAME ut_log-stats COMMAND test_ut_log-stats)
add_test (NAME ut_readMesh COMMAND test_ut_readMesh)
