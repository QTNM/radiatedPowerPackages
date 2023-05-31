#include "ElectronDynamics/PenningTraps.h"

rad::IdealPenningTrap::IdealPenningTrap(double BField, double v0, double rho0,
                                        double z0)
    : B(TVector3(0, 0, BField)), V0(v0), Rho0(rho0), Z0(z0) {}

TVector3 rad::IdealPenningTrap::evaluate_field_at_point(const TVector3 vec) {
  return B;
}

TVector3 rad::IdealPenningTrap::evaluate_e_field_at_point(TVector3 v) {
  double c{2 * V0 / (Z0 * Z0 + 0.5 * Rho0 * Rho0)};
  return c * TVector3(0.5 * v.X(), 0.5 * v.Y(), -v.Z());
}