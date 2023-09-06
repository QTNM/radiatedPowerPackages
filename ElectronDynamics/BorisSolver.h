/*
  BorisSolver.h

  Class containing an implementation of a Boris push method for solving electron
  dynamics
*/

#ifndef BORIS_SOLVER_H
#define BORIS_SOLVER_H

#include <tuple>

#include "BasicFunctions/BasicFunctions.h"
#include "BasicFunctions/Constants.h"
#include "ElectronDynamics/BaseField.h"
#include "TMath.h"
#include "TVector3.h"
#include "Waveguides/CircularCavity.h"

namespace rad {
class BorisSolver {
 private:
  double mass;
  double charge;
  double tau;
  BaseField* field = 0;
  CircularCavity* cav = 0;

  /// @brief Calculate omega as a function of position
  /// @param pos Position of charge
  TVector3 get_omega(const TVector3 pos);

  /// @brief Calculate radiation acceleration from position and velocity
  /// @param pos Position of the charge
  /// @param vec Velocity of the charge
  /// @return a 3-vector of the acceleration from the RR force
  TVector3 radiation_acceleration(const TVector3 pos, const TVector3 vel);

  /// @brief Returns the B field at the position
  /// @param pos Electron position
  /// @return Magnetic field vector in Tesla
  TVector3 calc_b_field(TVector3 pos);

  /// @brief Returns the E field at the position
  /// @param pos Electron position
  /// @return Magnetic field vector in volts/metre
  TVector3 calc_e_field(TVector3 pos);

 public:
  /// @brief Default constructor
  /// By default a uniform field is used
  BorisSolver();

  /// @brief Parametrised constructor
  /// @param field_v Pointer to a magnetic field calculator
  /// @param charge_v Particle charge. Default is electron charge
  /// @param mass_v Particle mass. Default is electron mass
  /// @param tau_v Energy loss. Default is zero
  /// @param cavity Pointer to resonant cavity (if using)
  BorisSolver(BaseField* field_v, const double charge_v = -TMath::Qe(),
              const double mass_v = ME, const double tau_v = 0.0,
              CircularCavity* cavity = 0);

  /// @brief Advances position and velocity vector by a set time
  /// @param time_step The time step over which to advance the charge's motion
  /// @param x0 Vector of the charge's starting position
  /// @param v0 Vector of the charge's starting velocity
  /// @return Tuple containing (1) the output position vector and (2) the
  /// output velocity vector
  std::tuple<TVector3, TVector3> advance_step(const double time_step,
                                              const TVector3 x0,
                                              const TVector3 v0);

  /// @brief acceleration due to B field and RR force
  /// @param pos Position of the charge
  /// @param vec Velocity of the charge
  /// @return a 3-vector of the acceleration
  TVector3 acc(const TVector3 pos, const TVector3 vel);
};
}  // namespace rad

#endif
