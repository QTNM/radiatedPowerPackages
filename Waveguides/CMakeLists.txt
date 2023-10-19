add_library(Waveguides IWaveguide.cxx CircularWaveguide.cxx RectangularWaveguide.cxx ICavity.cxx CircularCavity.cxx)
target_link_libraries(Waveguides PUBLIC BasicFunctions ${ROOT_LIBRARIES} ${Boost_MATH_LIBRARY})
target_include_directories(Waveguides PRIVATE ${Boost_INCLUDE_DIRS})