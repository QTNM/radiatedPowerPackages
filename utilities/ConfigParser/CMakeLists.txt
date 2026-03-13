add_library(ConfigParser
    ConfigParser.cxx
    FieldFactory.cxx
    DetectorFactory.cxx
    ElectronConfig.cxx
    ElectronGenerator.cxx
    SimulationConfig.cxx
)

target_link_libraries(ConfigParser
    PUBLIC  BasicCore ElectronDynamics Antennas Waveguides yaml-cpp ${ROOT_LIBRARIES}
)

target_include_directories(ConfigParser PUBLIC
    $<BUILD_INTERFACE:${CMAKE_SOURCE_DIR}>
    $<INSTALL_INTERFACE:include>
)
