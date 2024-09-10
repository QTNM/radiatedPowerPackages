#include "ElectronDynamics/ManyParticleSolver.h"

TVector3 rad::ManyParticleSolver::get_omega(Particle &p) {
  TVector3 BField{field->evaluate_field_at_point(p.GetPosition())};
  return calculate_omega(BField, p.GetCharge(), 0.0, p.GetMass());
}

TVector3 rad::ManyParticleSolver::radiation_acceleration(Particle &p) {
  TVector3 omega{get_omega(p)};
  double denom{1 + tau * tau * omega.Dot(omega)};
  double accX{0};
  double accY{0};
  double accZ{0};

  TVector3 vel{p.GetVelocity()};
  accX -= tau * (omega.Z() * omega.Z() + omega.Y() * omega.Y()) * vel.X();
  accX += tau * omega.X() * (omega.Z() * vel.Z() + omega.Y() * vel.Y());

  accY -= tau * (omega.Z() * omega.Z() + omega.X() * omega.X()) * vel.Y();
  accY += tau * omega.Y() * (omega.Z() * vel.Z() + omega.X() * vel.X());

  accZ -= tau * (omega.X() * omega.X() + omega.Y() * omega.Y()) * vel.Z();
  accZ += tau * omega.Z() * (omega.X() * vel.X() + omega.Y() * vel.Y());

  TVector3 acc{accX / denom, accY / denom, accZ / denom};
  return acc;
}

TVector3 rad::ManyParticleSolver::calc_b_field(Particle &p) {
  return field->evaluate_field_at_point(p.GetPosition());
}

TVector3 rad::ManyParticleSolver::calc_e_field(Particle &p) {
  return field->evaluate_e_field_at_point(p.GetPosition());
}

TVector3 rad::ManyParticleSolver::acc(Particle &p) {
  // Lorentz force
  TVector3 omega{get_omega(p)};
  TVector3 acc{p.GetVelocity().Cross(omega)};

  // Add Larmor terms
  acc += radiation_acceleration(p);

  // Acceleration from electric field
  acc += p.GetCharge() * calc_e_field(p);

  // Now calculate the contributions from the electric fields of all the other
  // particles
  for (auto &particle : particles) {
    if (particle.GetID() != p.GetID()) {
      ROOT::Math::XYZVector pField{
          CalcEField(p.GetPosition(), particle.GetPosition(),
                     particle.GetVelocity(), particle.GetAcceleration())};
      TVector3 pFieldv{pField.X(), pField.Y(), pField.Z()};
      acc += p.GetCharge() * pFieldv;
    }
  }

  return acc;
}

void rad::ManyParticleSolver::AdvanceParticles(double dt) {}