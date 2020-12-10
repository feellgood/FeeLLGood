PROJECT (FEELLGOOD)

find_package (Boost COMPONENTS system filesystem unit_test_framework REQUIRED)

include_directories (${CMAKE_SOURCE_DIR}
  ${Boost_INCLUDE_DIRS}
  )

set( BOOST_TEST_LOG_LEVEL "all" )

add_definitions (-DBOOST_TEST_DYN_LINK)

add_executable (test_ut_pt3D ut_pt3D.cpp)
add_executable (test_ut_node ut_node.cpp)

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
