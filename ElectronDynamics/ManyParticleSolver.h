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
  /// @param p The particle in question
  /// @return Omega vector
  TVector3 get_omega(Particle &p);

  /// @brief Calculate radiation acceleration from position and velocity
  /// @param p The particle in question
  /// @return A 3-vector of the acceleration from the RR force
  TVector3 radiation_acceleration(Particle &p);

  /// @brief Returns the B field at the position
  /// @param p Particle
  /// @return Magnetic field vector in Tesla
  TVector3 calc_b_field(Particle &p);

  /// @brief Returns the E field at the particle position
  /// @param p Particle
  /// @return Electric field vector in volts/metre
  TVector3 calc_e_field(Particle &p);

  /// @brief Calculate the acceleration of a particle under the various forces
  /// @param p Particle in question
  /// @return Acceleration vector in m/s^2
  TVector3 acc(Particle &p);

 public:
  /// @brief Parameterised constructor
  /// @param field_v Pointer to a magnetic field instance
  /// @param particles_v Vector of particles to be tracked
  ManyParticleSolver(BaseField *field_v, std::vector<Particle> particles_v)
      : field{field_v}, particles{particles_v} {}

  /// @brief Getter for particle states
  /// @return Vector of particle states
  std::vector<Particle> GetParticles() { return particles; }

  /// @brief Advance the particle states in time
  /// @param dt Time step in seconds
  void AdvanceParticles(double dt);
};
}  // namespace rad

#endif