/*
  Event.h

  Takes as its input one or more electrons in a magnetic field, propagates them, calculates their
  induced voltages on one or more antennas and performs some basic signal processing
*/

#ifndef EVENT_H
#define EVENT_H

#include "EventGeneration/ParticleState.h"
#include "EventGeneration/FieldStorer.h"
#include "EventGeneration/OutputVars.h"
#include "ElectronDynamics/BaseField.h"
#include "ElectronDynamics/BorisSolver.h"
#include "Antennas/IAntenna.h"
#include "Antennas/AntennaArray.h"
#include "SignalProcessing/LocalOscillator.h"
#include "SignalProcessing/NoiseFunc.h"

#include "TTree.h"

#include <vector>

namespace rad
{
  class Event
  {
  public:
    /// Parametrised constructor
    /// \param particles A vector of the various particle states one wishes to propagate
    /// \param antennas A vector of the antennas at which signals are generated
    /// \param field The magnetic field in which the particles are propagated
    /// \param simStepSize The time step size to use in the simulation (in seconds)
    /// \param simTime Total time for particles to be simulated for (in seconds)
    Event(std::vector<ParticleState> particles, std::vector<IAntenna *> antennas,
          BaseField *field, double simStepSize, double simTime);

    /// Parametrised constructor
    /// \param particles A vector of the various particle states one wishes to propagate
    /// \param field The magnetic field in which the particles are propagated
    /// \param simStepSize The time step size to use in the simulation (in seconds)
    /// \param simTime Total time for particles to be simulated for (in seconds)
    Event(std::vector<ParticleState> particles, BaseField *field,
          double simStepSize, double simTime);

    /// Parametrised constructor using an antenna array
    /// \param particles A vector of the various particle states one wishes to propagate
    /// \param arrays A vector of the antennas arrays for which signals are generated
    /// \param field The magnetic field in which the particles are propagated
    /// \param simStepSize The time step size to use in the simulation (in seconds)
    /// \param simTime Total time for particles to be simulated for (in seconds)
    Event(std::vector<ParticleState> particles, std::vector<AntennaArray> arrays,
          BaseField *field, double simStepSize, double simTime);

    /// Propagates the particles for the time specified at construction
    /// \param outputFile The output file directory. Leave as default to write no ouput
    /// \param vars A vector containing the variables to be written to file
    void PropagateParticles(const char *outputFile = NULL, std::vector<OutputVar> vars = {});

    /// Gives the clock time of the event
    /// \return The Event clock time (in seconds)
    double GetClockTime() { return clockTime; }

    /// Get the state for one of the particles in the Event
    /// \param particleIndex Index of particle we wish to return
    /// \return A copy of the selected ParticleState
    ParticleState GetParticle(int particleIndex);

    /// Gets the number of particles in a given event, active or inactive
    /// \return  Number of particles in event
    unsigned int GetNParticles();

    /// Adds relevant info to perform any kind of signal processing
    /// \param osc The local oscillator to use for the mixing
    /// \param noise The noise to add to the signal
    /// \param sampleRate The sample rate to use (in Hertz)
    void AddSigProcInfo(LocalOscillator osc, std::vector<GaussianNoise> noise,
                        double sampleRate);

  private:
    /// Vector of particles
    std::vector<ParticleState> particleList;

    /// Vector of solvers corresponding to the particles
    std::vector<BorisSolver> solverList;

    /// Vector of antennas for signals to be generated on
    /// Each element of the array represents a separate antenna
    std::vector<IAntenna *> antennaList;

    // We can also do this with arrays of antenna elements
    std::vector<AntennaArray> arrayList;

    double clockTime;             // Lab clock time in seconds
    double simulationStepSize;    // Simulation step size in seconds
    double maximumSimulationTime; // How long do we want to simulate for
    BaseField *magneticField;     // The magnetic field for the simulation

    std::vector<unsigned int> antennaArrayMap;

    // These arrays are only used for writing outputs
    int nAntennas; // Number of separate antennas
    int nArrays;   // Number of arrays/separate readout channels

    static const int nMaxAntennas{30};

    double particleTime;
    double pos[3];
    double vel[3];
    double acc[3];
    double bField[3];
    double antEField[nMaxAntennas][3];
    double antBField[nMaxAntennas][3];
    double antVoltage[nMaxAntennas];

    // Variables for the signal processed data
    double sampleTime; // The next sample time
    double VI[nMaxAntennas]; // In phase voltage component
    double VQ[nMaxAntennas]; // Quadrature voltage component

    bool computeSigProc;

    // Signal processing things
    LocalOscillator localOsc;             // Local oscillator
    std::vector<GaussianNoise> noiseFunc; // Noise sources
    double sRate;                         // Sample rate (in Hertz) 

    ////////////////////// Private member functions ////////////////////////////

    /// Checks if a particle has started relative to the Event clock time
    /// \param part The particle which we are checking if it has started
    /// \return True if the particle should be being propagated
    bool ParticleStartCheck(ParticleState part);

    /// Moves a given particle forward one step in time (according to the supplied time step)
    /// \param particleNumber The index of the particle in the event to be propagated
    void AdvanceParticleStep(int particleNumber);

    /// Calculates the light propagation time between a particle and an antenna point
    /// \param particle The particle in question
    /// \param antenna The selected antenna
    /// \return The light propagation time in seconds
    double GetPropagationTime(ParticleState particle, IAntenna *antenna);

    /// Creates an output tree with the desired variables in it
    /// \param vars Vector containing the desired output variables
    /// \return A TTree with the desired variables as branches
    TTree *CreateOutputTree(std::vector<OutputVar> vars);

    /// Creates an output tree for signal-processed outputs
    TTree *CreateSampleTree();

    /// Adds the data from the particle state to the TTree
    /// Specifically, data doesn't involve any kind of propagation time
    /// \param outputTree The tree to add the data to
    /// \param vars The variables to write to the TTree
    void AddLocalParticleData(TTree *outputTree, std::vector<OutputVar> vars);

    /// Calculates the electric field from a particle at a specific position
    /// \param particleIndex The index of the particle to calculate
    /// \param antennaIndex The index of the antenna point to calculate
    /// \return The electric field vector with units of V/m
    TVector3 GetEFieldAtAntenna(unsigned int particleIndex,
                                unsigned int antennaIndex);

    /// Calculates the electric field from a particle at an antenna
    /// \param particleIndex The index of the particle to c#alculate
    /// \param ant Pointer to the antenna at which to calculate the field
    /// \return The electric field vector with units of V/m
    TVector3 GetEFieldAtAntenna(unsigned int particleIndex, IAntenna *ant);

    /// Calculate the magnetic field from a particle at a specific position
    /// \param particleIndex The index of the particle to calculate
    /// \param antennaIndex The index of the antenna point to calculate
    /// \return The magnetic field vector with units of Tesla
    TVector3 GetBFieldAtAntenna(unsigned int particleIndex,
                                unsigned int antennaIndex);

    /// Calculates the magnetic field from a particle at an antenna
    /// \param particleIndex The index of the particle to c#alculate
    /// \param ant Pointer to the antenna at which to calculate the field
    /// \return The magnetic field vector with units of T
    TVector3 GetBFieldAtAntenna(unsigned int particleIndex, IAntenna *ant);

    /// Checks if an array of output vars contains a specific one
    /// \param vars Vector of output variables to be checked
    /// \param testVar Specific variable to be checked for
    /// \return True if vector contains specific output variable
    bool VectorContainsVar(std::vector<OutputVar> vars, OutputVar testVar);

    /// Sets 2D arrays to (defined) nonsense values
    void ResetArrays();
  };
}

#endif
