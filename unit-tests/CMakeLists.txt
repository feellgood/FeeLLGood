PROJECT (FEELLGOOD)

find_package (Boost COMPONENTS system filesystem unit_test_framework REQUIRED)

find_package(TBB REQUIRED)
if(TBB_FOUND)
    message(STATUS "TBB Found")
else()
    message(FATAL_ERROR "TBB NOT FOUND")
endif()

include_directories (${CMAKE_SOURCE_DIR}
  ${Boost_INCLUDE_DIRS}
  )

set( BOOST_TEST_LOG_LEVEL "all" )

add_definitions (-DBOOST_TEST_DYN_LINK)

add_executable (test_ut_pt3D ut_pt3D.cpp)
add_executable (test_ut_node ut_node.cpp)

SET(SOURCES ../facette.cpp ../tetra.cpp ut_element.cpp)
add_executable (test_ut_element ${SOURCES})

SET(SOURCES ../facette.cpp ut_facette.cpp)
add_executable (test_ut_facette ${SOURCES})

SET(SOURCES ../tetra.cpp ut_tetra.cpp)
add_executable (test_ut_tetra ${SOURCES})

SET(SOURCES ../tetra.cpp ut_anisotropy.cpp)
add_executable (test_ut_anisotropy ${SOURCES})

add_executable (test_ut_tiny ut_tiny.cpp)

add_executable (test_ut_time_int ut_time_int.cpp)

add_executable (test_ut_log-stats ut_log-stats.cpp)

target_link_libraries(test_ut_pt3D
  ${Boost_FILESYSTEM_LIBRARY}
  ${Boost_SYSTEM_LIBRARY}
  ${Boost_UNIT_TEST_FRAMEWORK_LIBRARY}
  )

target_link_libraries(test_ut_node
  ${Boost_FILESYSTEM_LIBRARY}
  ${Boost_SYSTEM_LIBRARY}
  ${Boost_UNIT_TEST_FRAMEWORK_LIBRARY}
  )

target_link_libraries(test_ut_element
  ${Boost_FILESYSTEM_LIBRARY}
  ${Boost_SYSTEM_LIBRARY}
  ${Boost_UNIT_TEST_FRAMEWORK_LIBRARY}
  TBB::tbb
  )

target_link_libraries(test_ut_facette
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
  
target_link_libraries(test_ut_tiny
  ${Boost_FILESYSTEM_LIBRARY}
  ${Boost_SYSTEM_LIBRARY}
  ${Boost_UNIT_TEST_FRAMEWORK_LIBRARY}
  )

target_link_libraries(test_ut_time_int
  ${Boost_FILESYSTEM_LIBRARY}
  ${Boost_SYSTEM_LIBRARY}
  ${Boost_UNIT_TEST_FRAMEWORK_LIBRARY}
  )

target_link_libraries(test_ut_log-stats
  ${Boost_UNIT_TEST_FRAMEWORK_LIBRARY}
  )
#enable_testing()

add_test (NAME ut_pt3D_arithmetic COMMAND test_ut_pt3D)
add_test (NAME ut_node COMMAND test_ut_node)
add_test (NAME ut_element COMMAND test_ut_element)
add_test (NAME ut_facette COMMAND test_ut_facette)
add_test (NAME ut_tetra COMMAND test_ut_tetra)
add_test (NAME ut_anisotropy COMMAND test_ut_anisotropy)
add_test (NAME ut_tiny COMMAND test_ut_tiny)
add_test (NAME ut_time_int COMMAND test_ut_time_int)
add_test (NAME ut_log-stats COMMAND test_ut_log-stats)
