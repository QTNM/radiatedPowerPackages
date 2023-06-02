add_library(SignalProcessing NoiseFunc.cxx LocalOscillator.cxx Signal.cxx InducedVoltage.cxx LockInAmplifier.cxx ScaledSignal.cxx SignalQuick.cxx)
target_link_libraries(SignalProcessing PUBLIC ${ROOT_LIBRARIES} ${ROOTFFTW_LIB} BasicFunctions FieldClasses Antennas)
