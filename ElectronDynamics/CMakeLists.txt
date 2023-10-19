add_library(ElectronDynamics BaseField.cxx QTNMFields.cxx BorisSolver.cxx TrajectoryGen.cxx ComsolFields.cxx PenningTraps.cxx)
target_link_libraries(ElectronDynamics PUBLIC BasicFunctions Waveguides ${ROOT_LIBRARIES} ${Boost_MATH_LIBRARY})
target_include_directories(ElectronDynamics PRIVATE ${Boost_INCLUDE_DIRS})