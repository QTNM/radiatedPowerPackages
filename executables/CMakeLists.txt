# Add any executables
add_executable(sineTest sineTest.cxx)
target_link_libraries(sineTest PUBLIC ${EXTRA_LIBS} ${ROOT_LIBRARIES} ${ROOTFFTW_LIB})

add_executable(normTest normTest.cxx)
target_link_libraries(normTest PUBLIC ${EXTRA_LIBS} ${ROOT_LIBRARIES} ${ROOTFFTW_LIB})

add_executable(noiseTests noiseTests.cxx)
target_link_libraries(noiseTests PUBLIC ${EXTRA_LIBS} ${ROOT_LIBRARIES} ${ROOTFFTW_LIB})

add_executable(signalSum signalSum.cxx)
target_link_libraries(signalSum PUBLIC ${EXTRA_LIBS} ${ROOT_LIBRARIES} ${ROOTFFTW_LIB})

add_executable(perpendicularAntenna perpendicularAntenna.cxx)
target_link_libraries(perpendicularAntenna PUBLIC ${EXTRA_LIBS} ${ROOT_LIBRARIES} ${ROOTFFTW_LIB})

add_executable(rayleighDistributionChecks rayleighDistributionChecks.cxx)
target_link_libraries(rayleighDistributionChecks PUBLIC BasicFunctions SignalProcessing Antennas ${ROOT_LIBRARIES})

add_executable(writeTrajectory writeTrajectory.cxx)
target_link_libraries(writeTrajectory PUBLIC ElectronDynamics ${ROOT_LIBRARIES})

add_executable(idealCaseSpectrogram idealCaseSpectrogram.cxx)
target_link_libraries(idealCaseSpectrogram PUBLIC SignalProcessing Antennas ${ROOT_LIBRARIES})

add_executable(singleBinPower singleBinPower.cxx)
target_link_libraries(singleBinPower PUBLIC SignalProcessing Antennas BasicFunctions ${ROOT_LIBRARIES} ${ROOTFFTW_LIB})

add_executable(crudeEnergyResolution crudeEnergyResolution.cxx)
target_link_libraries(crudeEnergyResolution PUBLIC ElectronDynamics SignalProcessing Antennas BasicFunctions ${ROOT_LIBRARIES})

add_executable(harmonicTrapDetectability harmonicTrapDetectability.cxx)
target_link_libraries(harmonicTrapDetectability PUBLIC Antennas SignalProcessing ElectronDynamics BasicFunctions ${ROOT_LIBRARIES})

add_executable(PhasedArrayVariableField PhasedArrayVariableField.cxx)
target_link_libraries(PhasedArrayVariableField PUBLIC Antennas SignalProcessing ElectronDynamics BasicFunctions ${ROOT_LIBRARIES})

add_executable(RFBandwidth RFBandwidth.cxx)
target_link_libraries(RFBandwidth PRIVATE ElectronDynamics BasicFunctions ${ROOT_LIBRARIES})

add_executable(SingleAntennaPower SingleAntennaPower.cxx)
target_link_libraries(SingleAntennaPower PRIVATE ElectronDynamics BasicFunctions Antennas SignalProcessing FieldClasses ${ROOT_LIBRARIES})

add_executable(RadialEnergyResolution RadialEnergyResolution.cxx)
target_link_libraries(RadialEnergyResolution PRIVATE ElectronDynamics BasicFunctions ${ROOT_LIBRARIES})

add_executable(RealisticFields RealisticFields.cxx)
target_link_libraries(RealisticFields PUBLIC Antennas SignalProcessing ElectronDynamics BasicFunctions ${ROOT_LIBRARIES})

add_executable(RFBandwidthCRESDA RFBandwidthCRESDA.cxx)
target_link_libraries(RFBandwidthCRESDA PRIVATE ElectronDynamics BasicFunctions ${ROOT_LIBRARIES})