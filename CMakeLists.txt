cmake_minimum_required(VERSION 3.18)
project(radiatedPower VERSION 0.1.0)

# Compiler/build settings
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# Locate the ROOT package
list(APPEND CMAKE_PREFIX_PATH $ENV{ROOTSYS})

find_package(ROOT REQUIRED COMPONENTS MathCore RIO Hist Tree Net)
include(${ROOT_USE_FILE})
include_directories(${CMAKE_SOURCE_DIR} ${ROOT_INCLUDE_DIRS})
add_definitions(${ROOT_CXX_FLAGS})

# Locate Boost libraries
find_package(Boost 1.73.0 REQUIRED ALL HINTS $ENV{BOOST_DIR})
include_directories(${Boost_INCLUDE_DIRS})

# Locate HDF5
find_package(HDF5 COMPONENTS C CXX HL REQUIRED)

# Locate FFTW3
find_package(FFTW3 REQUIRED)

# Set build locations for libraries and binaries
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)

# Add libraries
add_subdirectory(BasicFunctions)

add_subdirectory(FieldClasses)

add_subdirectory(SignalProcessing)

add_subdirectory(Antennas)

add_subdirectory(ElectronDynamics)

add_subdirectory(Waveguides)

add_subdirectory(Scattering)

add_subdirectory(executables)