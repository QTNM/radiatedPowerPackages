add_library(Waveguides IWaveguide.cxx CircularWaveguide.cxx)
target_link_libraries(Waveguides PUBLIC BasicFunctions ${ROOT_LIBRARIES} ${Boost_MATH_LIBRARY})
