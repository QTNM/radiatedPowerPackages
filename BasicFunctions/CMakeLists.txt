add_library(BasicFunctions BasicFunctions.cxx EMFunctions.cxx TritiumSpectrum.cxx ComplexVector3.cxx)
target_link_libraries(BasicFunctions PUBLIC BasicCore SignalUtils ${ROOT_LIBRARIES})
