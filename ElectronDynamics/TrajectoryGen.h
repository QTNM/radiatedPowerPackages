/*
  TrajectoryGen.h

  Class designed to allow easy generation of electron trajectories
  These trajectories should be able to be read by the other classes by default
*/

#ifndef TRAJECTORY_GEN_H
#define TRAJECTORY_GEN_H

#include "ElectronDynamics/BaseField.h"
#include "ElectronDynamics/BorisSolver.h"

#include "BasicFunctions/Constants.h"

#include "TString.h"
#include "TVector3.h"

namespace rad
{
  class ElectronTrajectoryGen
  {
  private:
    BorisSolver solver;
    TString outputFilePath;
    double stepSize;
    double startTime;
    double simulationTime;
    TVector3 initialPosition;
    TVector3 initialVelocity;

  public:
    /// Parametrised constructor
    /// \param outputFile The output root file path
    /// \param The magnetic field map to be used
    /// \param initPos The initial electron position
    /// \param initVel The initial electron velocity
    /// \param simStepSize The simulation step size to use in seconds
    /// \param simStepSize The time to simulate in seconds
    /// \param initialSimTime The initial time in the simulation. Default is 0
    /// \param tau Energy loss. Default is 0
    ElectronTrajectoryGen(TString outputFile, BaseField *field,
                          TVector3 initPos, TVector3 initVel,
                          double simStepSize, double simTime,
                          double initialSimTime = 0.0, double tau = 0.0); 

    /// Generates the trajectory with the specified parameters
    void GenerateTraj();
  };
}

#endif
