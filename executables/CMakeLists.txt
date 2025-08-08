# Add any executables
add_executable(sineTest sineTest.cxx)
target_link_libraries(sineTest PRIVATE ROOTUtils BasicCore ${ROOT_LIBRARIES})

add_executable(normTest normTest.cxx)
target_link_libraries(normTest PRIVATE BasicCore ROOTUtils FieldClasses ${ROOT_LIBRARIES})

add_executable(noiseTests noiseTests.cxx)
target_link_libraries(noiseTests PRIVATE Antennas BasicCore ROOTUtils SignalProcessing ${ROOT_LIBRARIES})

add_executable(signalSum signalSum.cxx)
target_link_libraries(signalSum PRIVATE Antennas SignalProcessing ${ROOT_LIBRARIES})

add_executable(perpendicularAntenna perpendicularAntenna.cxx)
target_link_libraries(perpendicularAntenna PRIVATE Antennas SignalProcessing ${ROOT_LIBRARIES})

add_executable(rayleighDistributionChecks rayleighDistributionChecks.cxx)
target_link_libraries(rayleighDistributionChecks PUBLIC SignalProcessing Antennas ${ROOT_LIBRARIES})

add_executable(writeTrajectory writeTrajectory.cxx)
target_link_libraries(writeTrajectory PUBLIC ElectronDynamics ${ROOT_LIBRARIES})

add_executable(singleBinPower singleBinPower.cxx)
target_link_libraries(singleBinPower PUBLIC SignalProcessing Antennas ${ROOT_LIBRARIES})

add_executable(crudeEnergyResolution crudeEnergyResolution.cxx)
target_link_libraries(crudeEnergyResolution PUBLIC ElectronDynamics SignalProcessing Antennas ${ROOT_LIBRARIES})

add_executable(harmonicTrapDetectability harmonicTrapDetectability.cxx)
target_link_libraries(harmonicTrapDetectability PUBLIC Antennas SignalProcessing ElectronDynamics ${ROOT_LIBRARIES})

add_executable(WR42Jackson WR42Jackson.cxx)
target_link_libraries(WR42Jackson PUBLIC ElectronDynamics Waveguides ${ROOT_LIBRARIES})

add_executable(PhasedArrayVariableField PhasedArrayVariableField.cxx)
target_link_libraries(PhasedArrayVariableField PUBLIC Antennas SignalProcessing ElectronDynamics ${ROOT_LIBRARIES})

add_executable(RFBandwidth RFBandwidth.cxx)
target_link_libraries(RFBandwidth PRIVATE ElectronDynamics ${ROOT_LIBRARIES})

add_executable(SingleAntennaPower SingleAntennaPower.cxx)
target_link_libraries(SingleAntennaPower PRIVATE ElectronDynamics Antennas SignalProcessing FieldClasses ${ROOT_LIBRARIES})

add_executable(RadialEnergyResolution RadialEnergyResolution.cxx)
target_link_libraries(RadialEnergyResolution PRIVATE ElectronDynamics BasicCore ROOTUtils ${ROOT_LIBRARIES})

add_executable(RealisticFields RealisticFields.cxx)
target_link_libraries(RealisticFields PUBLIC Antennas SignalProcessing ElectronDynamics ${ROOT_LIBRARIES})

add_executable(RFBandwidthCRESDA RFBandwidthCRESDA.cxx)
target_link_libraries(RFBandwidthCRESDA PRIVATE ElectronDynamics ${ROOT_LIBRARIES})

add_executable(SignalRequirements SignalRequirements.cxx)
target_link_libraries(SignalRequirements PRIVATE BasicCore ROOTUtils ${ROOT_LIBRARIES})

add_executable(FakeTrackGen FakeTrackGen.cxx)
target_link_libraries(FakeTrackGen PRIVATE BasicCore ROOTUtils ${ROOT_LIBRARIES})

add_executable(PenningTrapTraj PenningTrapTraj.cxx)
target_link_libraries(PenningTrapTraj PRIVATE ElectronDynamics ${ROOT_LIBRARIES})

add_executable(AxialFreqPT AxialFreqPT.cxx)
target_link_libraries(AxialFreqPT PRIVATE ElectronDynamics ${ROOT_LIBRARIES})

add_executable(StartPositionPT StartPositionPT.cxx)
target_link_libraries(StartPositionPT PRIVATE ElectronDynamics ${ROOT_LIBRARIES})

add_executable(NicePenningTraj NicePenningTraj.cxx)
target_link_libraries(NicePenningTraj PRIVATE ElectronDynamics ${ROOT_LIBRARIES})

add_executable(Project8CCA Project8CCA.cxx)
target_link_libraries(Project8CCA PRIVATE BasicCore ROOTUtils Waveguides ${ROOT_LIBRARIES})

add_executable(CRESDA0Cavity CRESDA0Cavity.cxx)
target_link_libraries(CRESDA0Cavity PRIVATE BasicCore ROOTUtils Waveguides ${ROOT_LIBRARIES})

add_executable(CavityFieldAmps CavityFieldAmps.cxx)
target_link_libraries(CavityFieldAmps PRIVATE Waveguides ElectronDynamics ${ROOT_LIBRARIES})

add_executable(ScatteringCheck ScatteringCheck.cxx)
target_link_libraries(ScatteringCheck PRIVATE ElectronDynamics Scattering ${ROOT_LIBRARIES})

add_executable(LongTimeGeneration LongTimeGeneration.cxx)
target_link_libraries(LongTimeGeneration PRIVATE ElectronDynamics Scattering ${ROOT_LIBRARIES})

add_executable(VariableEnergyLoss VariableEnergyLoss.cxx)
target_link_libraries(VariableEnergyLoss PRIVATE ElectronDynamics ${ROOT_LIBRARIES})

add_executable(HelmholtzCoils HelmholtzCoils.cxx)
target_link_libraries(HelmholtzCoils PRIVATE ElectronDynamics ${ROOT_LIBRARIES} ${Boost_MATH_LIBRARY})

add_executable(LeeWhitingCoils LeeWhitingCoils.cxx)
target_link_libraries(LeeWhitingCoils PRIVATE ElectronDynamics ${ROOT_LIBRARIES} ${Boost_MATH_LIBRARY})

add_executable(TrapCoilTuning TrapCoilTuning.cxx)
target_link_libraries(TrapCoilTuning PRIVATE ElectronDynamics ${ROOT_LIBRARIES})

add_executable(FourCoilEnergyRes FourCoilEnergyRes.cxx)
target_link_libraries(FourCoilEnergyRes PRIVATE ElectronDynamics SignalProcessing Waveguides ${ROOT_LIBRARIES})

add_executable(SignalFromWaveguide SignalFromWaveguide.cxx)
target_link_libraries(SignalFromWaveguide PRIVATE ElectronDynamics SignalProcessing Waveguides ${ROOT_LIBRARIES})

add_executable(P8HarmonicTrap P8HarmonicTrap.cxx)
target_link_libraries(P8HarmonicTrap PRIVATE ElectronDynamics SignalProcessing Waveguides ${ROOT_LIBRARIES})

add_executable(WaveguideSignalGen WaveguideSignalGen.cxx)
target_link_libraries(WaveguideSignalGen PRIVATE ElectronDynamics SignalProcessing Waveguides ${ROOT_LIBRARIES} ${HDF5_CXX_LIBRARIES})
target_include_directories(WaveguideSignalGen PRIVATE ${HDF5_INCLUDE_DIRS})

add_executable(P8WaveguidePower P8WaveguidePower.cxx)
target_link_libraries(P8WaveguidePower PRIVATE ElectronDynamics SignalProcessing Waveguides ${ROOT_LIBRARIES})

add_executable(EventRate EventRate.cxx)
target_link_libraries(EventRate PRIVATE BasicCore BasicFunctions ${ROOT_LIBRARIES})

add_executable(LargeVolumeCavityP8 LargeVolumeCavityP8.cxx)
target_link_libraries(LargeVolumeCavityP8 PRIVATE Waveguides ${ROOT_LIBRARIES})

add_executable(SyntheticSignals SyntheticSignals.cxx)
target_link_libraries(SyntheticSignals PRIVATE Antennas ElectronDynamics FieldClasses ${ROOT_LIBRARIES})

add_executable(TwoElectronSim TwoElectronSim.cxx)
target_link_libraries(TwoElectronSim PRIVATE ElectronDynamics ${ROOT_LIBRARIES})

add_executable(EventGenWithScattering EventGenWithScattering.cxx)
target_link_libraries(EventGenWithScattering PRIVATE ElectronDynamics Scattering SignalProcessing ${ROOT_LIBRARIES} ${HDF5_CXX_LIBRARIES})
target_include_directories(EventGenWithScattering PRIVATE ${HDF5_INCLUDE_DIRS})

add_executable(TemplateGeneration TemplateGeneration.cxx)
target_link_libraries(TemplateGeneration PRIVATE ElectronDynamics Scattering SignalProcessing ${ROOT_LIBRARIES} ${HDF5_CXX_LIBRARIES})
target_include_directories(TemplateGeneration PRIVATE ${HDF5_INCLUDE_DIRS})

add_executable(GenerateAtPitchAngle GenerateAtPitchAngle.cxx)
target_link_libraries(GenerateAtPitchAngle PRIVATE ElectronDynamics SignalProcessing Waveguides ${ROOT_LIBRARIES} ${HDF5_CXX_LIBRARIES})
target_include_directories(GenerateAtPitchAngle PRIVATE ${HDF5_INCLUDE_DIRS})

add_executable(TestWgTheory TestWgTheory.cxx)
target_link_libraries(TestWgTheory PRIVATE ElectronDynamics SignalProcessing ${ROOT_LIBRARIES} ${HDF5_CXX_LIBRARIES})
target_include_directories(TestWgTheory PRIVATE ${HDF5_INCLUDE_DIRS})

add_executable(MLSampleGeneration MLSampleGeneration.cxx)
target_link_libraries(MLSampleGeneration PRIVATE ElectronDynamics SignalProcessing Waveguides ${HDF5_CXX_LIBRARIES})
target_include_directories(MLSampleGeneration PRIVATE ${HDF5_INCLUDE_DIRS})

add_executable(SecondaryElectronProduction SecondaryElectronProduction.cxx)
target_link_libraries(SecondaryElectronProduction PRIVATE ElectronDynamics Scattering ${ROOT_LIBRARIES})
