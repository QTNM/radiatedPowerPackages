add_library(BasicFunctions BasicFunctions.cxx EMFunctions.cxx TritiumSpectrum.cxx ButterworthFilter.cxx)
target_link_libraries(BasicFunctions PUBLIC ${ROOT_LIBRARIES} ${ROOTFFTW_LIB})
