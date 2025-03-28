add_executable(testXSecs testXSecs.cxx)
target_link_libraries(testXSecs PRIVATE BasicFunctions Scattering ${ROOT_LIBRARIES})