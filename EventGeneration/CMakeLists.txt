add_library(EventGeneration ITritiumSpectrum.cxx BasicTritiumSpectrum.cxx ParticleState.cxx Event.cxx FieldStorer.cxx)
target_link_libraries(EventGeneration PUBLIC ${ROOT_LIBRARIES} BasicFunctions ElectronDynamics Antennas SignalProcessing)
