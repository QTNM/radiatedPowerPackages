add_library(ElectronDynamics 
  BaseField.cxx 
  QTNMFields.cxx 
  BorisSolver.cxx 
  TrajectoryGen.cxx 
  ComsolFields.cxx 
  PenningTraps.cxx
)

target_link_libraries(ElectronDynamics PUBLIC 
  BasicCore 
  Waveguides 
  ${ROOT_LIBRARIES} 
  ${Boost_MATH_LIBRARY}
)
