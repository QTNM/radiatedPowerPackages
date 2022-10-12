/*
  ParticlePoint.h

  Represents the instantaneous state of a particle at stored time
*/

#ifndef PARTICLE_STATE_H
#define PARTICLE_STATE_H

#include "TVector3.h"

namespace rad
{
  class ParticleState {
  public:
    /// Parametrised constructor
    /// \param tStart Particle start time in seconds
    /// \param mass Particle mass in kg
    /// \param charge Particle charge in coulombs
    /// \param initialPos Vector of particle initial position (in metres)
    /// \param initialVel Vector of particle initial velocity (in m/s)
    ParticleState(double tStart, double mass, double charge, TVector3 initialPos, TVector3 initialVel);

    /// Access the position vector at the current time
    /// \Returns The position vector at the current time
    TVector3 GetPositionVector() { return positionVector; }

    /// Access the velocity vector at the current time
    /// \Returns The velocity vector at the current time
    TVector3 GetVelocityVector() { return velocityVector; }

    /// Access the particle mass
    /// \Returns The particle mass (in kilograms) 
    double GetParticleMass() { return particleMass; }

    /// Access the particle charge
    /// \Returns The particle charge (in coulombs)
    double GetParticleCharge() { return particleCharge; }

    /// Access the current particle time
    /// \Returns The current simulation time of the particle (in seconds)
    double GetCurrentTime() { return currentTime; }

    /// Gets the relativistic energy of the particle
    /// \returns The energy of the particle in Joules
    double GetE();
    
    /// Gets the relativistic kinetic energy of the particle
    /// \returns The kinetic energy of the particle in Joules
    double GetKE();

  private:
    double startTime;        // Time in seconds for the particle to start propagating
    double currentTime;      // Time in seconds at which the particle has this state
    
    double particleMass;     // Particle mass (in kilograms)
    double particleCharge;   // Particle charge (in coulombs)
    
    TVector3 positionVector; // Position vector of the particle at the current time
    TVector3 velocityVector; // Velocity vector of the particle at the current time

    friend class Event;
  };
}

#endif
