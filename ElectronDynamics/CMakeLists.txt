add_library(ElectronDynamics BaseField.cxx QTNMFields.cxx BorisSolver.cxx TrajectoryGen.cxx ComsolFields.cxx PenningTraps.cxx)
target_link_libraries(ElectronDynamics PUBLIC BasicFunctions ${ROOT_LIBRARIES} ${Boost_MATH_LIBRARY})
