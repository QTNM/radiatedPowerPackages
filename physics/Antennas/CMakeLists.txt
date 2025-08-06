add_library(Antennas IAntenna.cxx HertzianDipole.cxx HalfWaveDipole.cxx PatchAntenna.cxx IsotropicAntenna.cxx)
target_link_libraries(Antennas PUBLIC BasicFunctions ${ROOT_LIBRARIES})
