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
find_package(HDF5 COMPONENTS C CXX HL REQUIRED HINTS $ENV{HDF5_INSTALL})
link_directories(${HDF5_LIBRARY_DIRS})
include_directories(${HDF5_INCLUDE_DIRS})

# Add the FFTW ROOT wrapper library
# Need to check if the helper variable is set
if(DEFINED ENV{FFTW_WRAPPER_INSTALL})
	message("Detected sugggested location for ROOT FFTW wrapper")
	set(FFTW_WRAPPER_LIB_DIR $ENV{FFTW_WRAPPER_INSTALL}/lib)
	set(FFTW_WRAPPER_INC_DIR $ENV{FFTW_WRAPPER_INSTALL}/include)
	find_library(ROOTFFTW_LIB RootFftwWrapper REQUIRED HINTS ${FFTW_WRAPPER_LIB_DIR})
	include_directories(${FFTW_WRAPPER_INC_DIR})
else()
	message("Looking for ROOT FFTW wrapper in default location")
	find_library(ROOTFFTW_LIB RootFftwWrapper REQUIRED)
endif()
message("ROOTFFTW_LIB is set to ${ROOTFFTW_LIB}")

# Set build locations for libraries and binaries
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)

# Add libraries
add_subdirectory(BasicFunctions)
list(APPEND EXTRA_LIBS BasicFunctions)
list(APPEND EXTRA_INCLUDES "${PROJECT_SOURCE_DIR}/BasicFunctions")

add_subdirectory(FieldClasses)
list(APPEND EXTRA_LIBS FieldClasses)
list(APPEND EXTRA_INCLUDES "${PROJECT_SOURCE_DIR}/FieldClasses")

add_subdirectory(SignalProcessing)
list(APPEND EXTRA_LIBS SignalProcessing)
list(APPEND EXTRA_INCLUDES "${PROJECT_SOURCE_DIR}/SignalProcessing")

add_subdirectory(Antennas)
list(APPEND EXTRA_LIBS Antennas)
list(APPEND EXTRA_INCLUDES "${PROJECT_SOURCE_DIR}/Antennas")

add_subdirectory(ElectronDynamics)

add_subdirectory(Scattering)

add_subdirectory(executables)

