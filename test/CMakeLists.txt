add_executable(testXSecs testXSecs.cxx)
target_link_libraries(testXSecs PRIVATE ROOTUtils BasicCore Scattering ${ROOT_LIBRARIES})

add_executable(testWaveguide testWaveguide.cxx)
target_link_libraries(testWaveguide PRIVATE  ElectronDynamics SignalProcessing Waveguides ${HDF5_CXX_LIBRARIES} ${ROOT_LIBRARIES})
target_include_directories(testWaveguide PRIVATE ${HDF5_INCLUDE_DIRS})