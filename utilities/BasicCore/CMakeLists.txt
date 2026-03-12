# BasicCore - Header-only library with no external dependencies
add_library(BasicCore INTERFACE)

# Set include directories for interface library
target_include_directories(BasicCore INTERFACE 
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/..>
    $<INSTALL_INTERFACE:include>
)

# No linking required - pure header-only library