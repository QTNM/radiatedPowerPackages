add_library(Waveguides IWaveguide.cxx CircularWaveguide.cxx RectangularWaveguide.cxx WaveguideMode.cxx ICavity.cxx CircularCavity.cxx Probe.cxx)
target_link_libraries(Waveguides PUBLIC BasicFunctions ${ROOT_LIBRARIES} ${Boost_MATH_LIBRARY})
