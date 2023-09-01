add_library(SignalProcessing NoiseFunc.cxx LocalOscillator.cxx Signal.cxx InducedVoltage.cxx)
target_link_libraries(SignalProcessing PUBLIC ${ROOT_LIBRARIES} ${ROOTFFTW_LIB} BasicFunctions FieldClasses Antennas)
