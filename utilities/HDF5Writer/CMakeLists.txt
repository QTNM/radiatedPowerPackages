add_library(HDF5Writer
    HDF5Writer.cxx
)

target_link_libraries(HDF5Writer
    PUBLIC ${HDF5_CXX_LIBRARIES} ${ROOT_LIBRARIES}
)

target_include_directories(HDF5Writer PUBLIC
    $<BUILD_INTERFACE:${CMAKE_SOURCE_DIR}>
    $<INSTALL_INTERFACE:include>
    ${HDF5_INCLUDE_DIRS}
)
