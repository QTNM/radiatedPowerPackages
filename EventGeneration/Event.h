/*
  Event.h

  Takes as its input one or more electrons in a magnetic field, propagates them, calculates their
  induced voltages on one or more antennas and performs some basic signal processing
*/

#ifndef EVENT_H
#define EVENT_H

#include "EventGeneration/ParticleState.h"
#include "ElectronDynamics/BaseField.h"
#include "ElectronDynamics/BorisSolver.h"

#include <vector>

namespace rad
{
  class Event {
  public:

    /// Parametrised constructor
    /// \param particles A vector of the various particle states one wishes to propagate
    /// \param field The magnetic field in which the particles are propagated
    /// \param simStepSize The time step size to use in the simulation (in seconds)
    /// \param simTime Total time for particles to be simulated for (in seconds)
    Event(std::vector<ParticleState> particles, BaseField* field, double simStepSize, double simTime);

    /// Propagates the particles for the time specified at construction
    void PropagateParticles();
    
  private:
    /// Vector of particles
    std::vector<ParticleState> particleList;

    /// Vector of solvers corresponding to the particles
    std::vector<BorisSolver> solverList;
    
    double clockTime;
    double simulationStepSize;
    double maximumSimulationTime;
    BaseField* magneticField;

    /// Checks if a particle has started relative to the Event clock time
    /// \param part The particle which we are checking if it has started
    /// \Returns True if the particle should be being propagated
    bool ParticleStartCheck(ParticleState part);

    /// Moves a given particle forward one step in time (according to the supplied time step)
    /// \param particleNumber The index of the particle in the event to be propagated
    void AdvanceParticleStep(int particleNumber);
  };
}

#endif
