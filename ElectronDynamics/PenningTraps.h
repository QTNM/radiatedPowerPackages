/*
  PenningTraps.h

  Derived field classes for Penning traps

  S. Jones 04-04-23
*/

#ifndef PENNING_TRAPS_H
#define PENNING_TRAPS_H

#include "ElectronDynamics/BaseField.h"
#include "TVector3.h"

namespace rad {
class IdealPenningTrap : public BaseField {
 public:
  /// @brief Parametrised constructor
  /// @param BField Magnetic field strength [Tesla]
  /// @param v0 Voltage between ring electron and endcaps [Volts]
  /// @param rho0 Distance between trap centre and ring electrode [m]
  /// @param z0 Distance between trap centre and endcaps [m]
  IdealPenningTrap(double BField, double v0, double rho0, double z0);

  /// @brief Calculate B field at a point
  /// @param vec Position vector
  /// @return B field vector in tesla
  TVector3 evaluate_field_at_point(const TVector3 vec) override;

  /// @brief Calculate E field at point
  /// @param v Position vector
  /// @return E field vector in volts / metre
  TVector3 evaluate_e_field_at_point(TVector3 v) override;

 private:
  TVector3 B;
  double V0;
  double Rho0;
  double Z0;
};
}  // namespace rad

#endif