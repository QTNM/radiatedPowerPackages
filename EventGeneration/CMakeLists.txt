add_library(EventGeneration ITritiumSpectrum.cxx BasicTritiumSpectrum.cxx)
target_link_libraries(EventGeneration PUBLIC ${ROOT_LIBRARIES} BasicFunctions)
