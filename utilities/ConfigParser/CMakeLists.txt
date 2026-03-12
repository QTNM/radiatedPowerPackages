add_library(ConfigParser
    ConfigParser.cxx
    FieldFactory.cxx
    ElectronConfig.cxx
    SimulationConfig.cxx
)

target_link_libraries(ConfigParser
    PUBLIC  BasicCore ElectronDynamics yaml-cpp::yaml-cpp ${ROOT_LIBRARIES}
)

target_include_directories(ConfigParser PUBLIC
    $<BUILD_INTERFACE:${CMAKE_SOURCE_DIR}>
    $<INSTALL_INTERFACE:include>
)
