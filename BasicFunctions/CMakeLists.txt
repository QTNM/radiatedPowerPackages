add_library(BasicFunctions BasicFunctions.cxx EMFunctions.cxx TritiumSpectrum.cxx ButterworthFilter.cxx ComplexVector3.cxx FFTWComplex.cxx FourierTransforms.cxx)
target_link_libraries(BasicFunctions PUBLIC BasicCore ${ROOT_LIBRARIES} ${FFTW3_LIBRARIES})
