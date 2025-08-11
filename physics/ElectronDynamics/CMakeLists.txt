add_library(ElectronDynamics 
  BaseField.cxx 
  QTNMFields.cxx 
  BorisSolver.cxx 
  TrajectoryGen.cxx 
  ComsolFields.cxx 
  PenningTraps.cxx
)

target_link_libraries(ElectronDynamics 
  PUBLIC 
  BasicCore 
  ${ROOT_LIBRARIES} 
  PRIVATE
  Waveguides 
  ${Boost_MATH_LIBRARY}
)
