# Utility Libraries
# These libraries provide reusable functionality across the project

# Header-only constants and basic math utilities
add_subdirectory(BasicCore)

# FFTW and Boost dependent signal processing utilities  
add_subdirectory(SignalUtils)

# ROOT framework dependent utilities
add_subdirectory(ROOTUtils)

# Legacy consolidated utility functions
add_subdirectory(BasicFunctions)

# YAML configuration parsing utilities
add_subdirectory(ConfigParser)