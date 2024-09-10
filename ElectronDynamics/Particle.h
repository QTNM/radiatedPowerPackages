/*
  Particle.h

  Class representing the state of a particle at a given time
*/

#ifndef PARTICLE_H
#define PARTICLE_H

#include "BasicFunctions/Constants.h"
#include "TMath.h"
#include "TVector3.h"

namespace rad {
class Particle {
 private:
  unsigned int id;
  static unsigned int counter;

  TVector3 position;
  TVector3 velocity;
  TVector3 acceleration;
  double mass;
  double charge;

 public:
  /// @brief Parameterised constructor
  /// @param pos Initial position vector of the particle in metres
  /// @param vel Initial velocity vector of the particle in m/s
  /// @param acc Initial acceleration vector of the particle in m/s^2
  /// @param m Mass of the particle in kg
  /// @param q Charge of the particle in C
  Particle(TVector3 pos, TVector3 vel, TVector3 acc, double m = ME,
           double q = -TMath::Qe());

  unsigned int GetID() { return id; }

  /// @brief Getter for particle position
  /// @return Particle position vector in metres
  TVector3 GetPosition() { return position; }

  /// @brief Getter for particle velocity
  /// @return Particle velocity vector in m/s
  TVector3 GetVelocity() { return velocity; }

  /// @brief Getter for particle acceleration
  /// @return Particle acceleration vector in m/s^2
  TVector3 GetAcceleration() { return acceleration; }

  /// @brief Getter for particle mass
  /// @return Particle mass in kg
  double GetMass() { return mass; }

  /// @brief Getter for particle charge
  /// @return Particle charge in C
  double GetCharge() { return charge; }

  /// @brief Function to update the state of the particle
  /// @param newpos New position vector of the particle in metres
  /// @param newvel New velocity vector of the particle in m/s
  /// @param newacc New acceleration vector of the particle in m/s^2
  void UpdateState(TVector3 newpos, TVector3 newvel, TVector3 newacc);
};
}  // namespace rad

#endif