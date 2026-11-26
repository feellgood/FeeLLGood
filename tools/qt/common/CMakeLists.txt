# Code shared by the feeLLGood tools.
add_library(fgcommon STATIC
    gmshMesh.cpp
    jsMath.cpp
    pyRepr.cpp
    textTable.cpp)
target_include_directories(fgcommon PUBLIC ${CMAKE_CURRENT_SOURCE_DIR} ${GMSH_INC})
target_link_libraries(fgcommon PUBLIC Qt5::Core Qt5::Qml ${GMSH_LIB})

# HDF5 helpers of fg-sol2hdf5.
if(HDF5_FOUND)
    add_library(fgh5 STATIC h5Utils.cpp)
    target_include_directories(fgh5 PUBLIC ${HDF5_INCLUDE_DIRS})
    target_compile_definitions(fgh5 PUBLIC ${HDF5_DEFINITIONS})
    target_link_libraries(fgh5 PUBLIC fgcommon ${HDF5_CXX_LIBRARIES} ${HDF5_LIBRARIES})
endif()
