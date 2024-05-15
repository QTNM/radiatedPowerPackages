add_library(SignalProcessing NoiseFunc.cxx LocalOscillator.cxx Signal.cxx InducedVoltage.cxx)
target_link_libraries(SignalProcessing PUBLIC ${ROOT_LIBRARIES} BasicFunctions FieldClasses Antennas Waveguides)
