/*
  BorisSolver.cxx

  Function implementations for the energy conserving Boris solver
*/

#include "ElectronDynamics/BorisSolver.h"

#include <tuple>

#include "ElectronDynamics/BaseField.h"
#include "ElectronDynamics/QTNMFields.h"
#include "TMath.h"
#include "TVector3.h"

rad::BorisSolver::BorisSolver()
    : mass(ME), charge(-TMath::Qe()), tau(0), field(new UniformField(1.0)) {}

rad::BorisSolver::BorisSolver(BaseField* field_v, const double charge_v,
                              const double mass_v, const double tau_v)
    : charge(charge_v), mass(mass_v), field(field_v), tau(tau_v) {}

TVector3 rad::BorisSolver::get_omega(const TVector3 pos) {
  TVector3 BField = field->evaluate_field_at_point(pos);
  return calculate_omega(BField, charge, 0.0, mass);
}

TVector3 rad::BorisSolver::radiation_acceleration(const TVector3 pos,
                                                  const TVector3 vel) {
  TVector3 omega = get_omega(pos);
  double denom = 1 + tau * tau * omega.Dot(omega);
  double accX = 0;
  double accY = 0;
  double accZ = 0;

  accX -= tau * (omega.Z() * omega.Z() + omega.Y() * omega.Y()) * vel.X();
  accX += tau * omega.X() * (omega.Z() * vel.Z() + omega.Y() * vel.Y());

  accY -= tau * (omega.Z() * omega.Z() + omega.X() * omega.X()) * vel.Y();
  accY += tau * omega.Y() * (omega.Z() * vel.Z() + omega.X() * vel.X());

  accZ -= tau * (omega.X() * omega.X() + omega.Y() * omega.Y()) * vel.Z();
  accZ += tau * omega.Z() * (omega.X() * vel.X() + omega.Y() * vel.Y());

  TVector3 acc(accX / denom, accY / denom, accZ / denom);
  return acc;
}

TVector3 rad::BorisSolver::acc(const TVector3 pos, const TVector3 vel) {
  TVector3 omega = get_omega(pos);

  // Lorentz force
  TVector3 acc = vel.Cross(omega);

  // Add Larmor terms
  acc += radiation_acceleration(pos, vel);

  return acc;
}

std::tuple<TVector3, TVector3> rad::BorisSolver::advance_step(
    const double time_step, const TVector3 x0, const TVector3 v0) {
  double gamma_n = 1 / sqrt(1 - v0.Dot(v0) / pow(TMath::C(), 2));
  TVector3 u_n = v0 * gamma_n;
  TVector3 x_n = x0;
  TVector3 v_n = v0;

  // Half position step
  TVector3 x_nplushalf = x_n + v_n * (time_step / 2.0);

  // Do the first half of the RR acceleration
  TVector3 u_minus =
      u_n + (time_step / 2.0) * radiation_acceleration(x_nplushalf, u_n);
  double gamma_minus = sqrt(1.0 + u_n.Dot(u_n) / pow(TMath::C(), 2));

  // Rotation step
  TVector3 B_nplushalf = calc_b_field(x_nplushalf);

  double theta = charge * time_step / (mass * gamma_minus) * B_nplushalf.Mag();
  TVector3 t = TMath::Tan(theta / 2.0) * B_nplushalf.Unit();
  TVector3 u_prime = u_minus + u_minus.Cross(t);
  TVector3 u_plus = u_minus + 2.0 / (1.0 + t.Dot(t)) * u_prime.Cross(t);

  // Second half of the RR acceleration
  TVector3 u_nplus1 =
      u_plus + (time_step / 2.0) * radiation_acceleration(x_nplushalf, u_plus);

  // Now update position
  double gamma_nplus1 = sqrt(1 + pow(u_nplus1.Mag() / TMath::C(), 2));
  TVector3 v_nplus1 = u_nplus1 * (1 / gamma_nplus1);
  TVector3 x_nplus1 = x_nplushalf + v_nplus1 * (time_step / 2.0);

  return std::make_tuple(x_nplus1, v_nplus1);
}

TVector3 rad::BorisSolver::calc_b_field(const TVector3 pos) {
  return (field->evaluate_field_at_point(pos));
}
