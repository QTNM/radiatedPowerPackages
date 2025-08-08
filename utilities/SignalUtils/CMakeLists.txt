# SignalUtils - Signal processing utilities with FFTW and Boost dependencies
add_library(SignalUtils 
    FFTWComplex.cxx 
    FourierTransforms.cxx 
    BandpassFilter.cxx 
    ButterworthFilter.cxx
)

# Link to required external libraries
target_link_libraries(SignalUtils PUBLIC ${FFTW3_LIBRARIES} ${Boost_MATH_LIBRARY})

# Include directories
target_include_directories(SignalUtils PUBLIC 
    $<BUILD_INTERFACE:${CMAKE_SOURCE_DIR}>
    $<INSTALL_INTERFACE:include>
)