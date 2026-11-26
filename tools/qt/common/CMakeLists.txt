# Code shared by the feeLLGood tools.
add_library(fgcommon STATIC
    gmshMesh.cpp
    pyRepr.cpp
    textTable.cpp)
target_include_directories(fgcommon PUBLIC ${CMAKE_CURRENT_SOURCE_DIR} ${GMSH_INC})
target_link_libraries(fgcommon PUBLIC Qt5::Core ${GMSH_LIB})
