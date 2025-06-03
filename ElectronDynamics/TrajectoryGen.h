/*
  TrajectoryGen.h

  Class designed to allow easy generation of electron trajectories
  These trajectories should be able to be read by the other classes by default
*/

#ifndef TRAJECTORY_GEN_H
#define TRAJECTORY_GEN_H

#include "BasicFunctions/Constants.h"
#include "ElectronDynamics/BaseField.h"
#include "ElectronDynamics/BorisSolver.h"
#include "TString.h"
#include "TVector3.h"
#include "Waveguides/CircularCavity.h"

namespace rad {
class ElectronTrajectoryGen {
 private:
  BorisSolver solver;

 public:
  /// @brief Parametrised constructor
  /// @param outputFile The output root file path
  /// @param The magnetic field map to be used
  /// @param initPos The initial electron position
  /// @param initVel The initial electron velocity
  /// @param simStepSize The simulation step size to use in seconds
  /// @param simTime The time to simulate in seconds
  /// @param energyLoss Boolean for choosing if energy loss is on
  /// @param initialSimTime The initial time in the simulation. Default is 0
  ElectronTrajectoryGen(TString outputFile, BaseField* field, TVector3 initPos,
                        TVector3 initVel, double simStepSize, double simTime,
                        bool energyLoss = true, CircularCavity* cav = 0,
                        double initialSimTime = 0.0);
};
}  // namespace rad

#endif
