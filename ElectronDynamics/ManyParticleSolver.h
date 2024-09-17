/*
  ManyParticleSolver.h

  Class using a Boris push method for solving the dynamics of multiple electrons
*/

#ifndef MANY_PARTICLE_SOLVER_H
#define MANY_PARTICLE_SOLVER_H

#include <vector>

#include "BasicFunctions/BasicFunctions.h"
#include "BasicFunctions/Constants.h"
#include "BasicFunctions/EMFunctions.h"
#include "ElectronDynamics/BaseField.h"
#include "ElectronDynamics/Particle.h"
#include "TVector3.h"

namespace rad {
class ManyParticleSolver {
 private:
  BaseField *field = 0;
  std::vector<Particle> particles;

  const double tau = 2 * R_E / (3 * TMath::C());

  /// @brief Calculate omega as a function of position
  /// @param p The array index particle in question
  /// @return Omega vector
  TVector3 get_omega(TVector3 pos);

  /// @brief Calculate omega as a function of position
  /// @param p The array index particle in question
  /// @return Omega vector
  TVector3 get_omega(size_t pID);

  /// @brief Calculate radiation acceleration from position and velocity
  /// @param pID The particle in question
  /// @return A 3-vector of the acceleration from the RR force
  TVector3 radiation_acceleration(size_t pID);

  /// @brief Calculate radiation acceleration from position and velocity
  /// @param pos Position of the charge
  /// @param vec Velocity of the charge
  /// @return a 3-vector of the acceleration from the RR force
  TVector3 radiation_acceleration(TVector3 pos, TVector3 vel);

  /// @brief Returns the B field at the position
  /// @param pID Particle
  /// @return Magnetic field vector in Tesla
  TVector3 calc_b_field(size_t pID);

  /// @brief Returns the E field at the particle position
  /// @param parts Vector of all the particles
  /// @param pID ID of particle to calculate field at
  /// @return Electric field vector in volts/metre
  TVector3 calc_e_field(std::vector<Particle> &parts, size_t pID);

  /// @brief Calculate the acceleration of a particle under the various forces
  /// @param parts Vector of particles we want to calculate the field from
  /// @param pID ID of Particle in question
  /// @return Acceleration vector in m/s^2
  TVector3 acc(std::vector<Particle> &parts, size_t pID);

 public:
  /// @brief Parameterised constructor
  /// @param field_v Pointer to a magnetic field instance
  /// @param particles_v Vector of particles to be tracked
  ManyParticleSolver(BaseField *field_v, std::vector<Particle> particles_v);

  /// @brief Getter for particle states
  /// @return Vector of particle states
  std::vector<Particle> GetParticles() { return particles; }

  /// @brief Getter for a single particle state
  /// @param pID ID of Particle in question
  /// @return Particle state
  Particle GetParticle(size_t pID) { return particles.at(pID); }

  /// @brief Advance the particle states in time
  /// @param dt Time step in seconds
  void AdvanceParticles(double dt);
};
}  // namespace rad

#endif