#include "ElectronDynamics/ManyParticleSolver.h"

rad::ManyParticleSolver::ManyParticleSolver(BaseField *field_v,
                                            std::vector<Particle> particles_v)
    : field{field_v}, particles{particles_v} {
  // Calculate and set each particle's accleration
  for (size_t iPart{0}; iPart < particles.size(); iPart++) {
    particles.at(iPart).SetAcceleration(acc(particles, iPart));
  }
}

TVector3 rad::ManyParticleSolver::get_omega(size_t pID) {
  TVector3 BField{
      field->evaluate_field_at_point(particles.at(pID).GetPosition())};
  return calculate_omega(BField, particles.at(pID).GetCharge(), 0.0,
                         particles.at(pID).GetMass());
}

TVector3 rad::ManyParticleSolver::get_omega(TVector3 pos) {
  TVector3 BField{field->evaluate_field_at_point(pos)};
  return calculate_omega(BField, -TMath::Qe(), 0.0, ME);
}

TVector3 rad::ManyParticleSolver::radiation_acceleration(size_t pID) {
  TVector3 omega{get_omega(pID)};
  double denom{1 + tau * tau * omega.Dot(omega)};
  double accX{0};
  double accY{0};
  double accZ{0};

  TVector3 vel{particles.at(pID).GetVelocity()};
  accX -= tau * (omega.Z() * omega.Z() + omega.Y() * omega.Y()) * vel.X();
  accX += tau * omega.X() * (omega.Z() * vel.Z() + omega.Y() * vel.Y());

  accY -= tau * (omega.Z() * omega.Z() + omega.X() * omega.X()) * vel.Y();
  accY += tau * omega.Y() * (omega.Z() * vel.Z() + omega.X() * vel.X());

  accZ -= tau * (omega.X() * omega.X() + omega.Y() * omega.Y()) * vel.Z();
  accZ += tau * omega.Z() * (omega.X() * vel.X() + omega.Y() * vel.Y());

  TVector3 acc{accX / denom, accY / denom, accZ / denom};
  return acc;
}

TVector3 rad::ManyParticleSolver::radiation_acceleration(TVector3 pos,
                                                         TVector3 vel) {
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

TVector3 rad::ManyParticleSolver::calc_b_field(size_t pID) {
  return field->evaluate_field_at_point(particles.at(pID).GetPosition());
}

TVector3 rad::ManyParticleSolver::calc_e_field(std::vector<Particle> &parts,
                                               size_t pID) {
  TVector3 staticField{
      field->evaluate_e_field_at_point(particles.at(pID).GetPosition())};
  TVector3 particleField(0, 0, 0);
  for (size_t iPart{0}; iPart < parts.size(); iPart++) {
    if (iPart != pID) {
      ROOT::Math::XYZVector pField{CalcEField(
          parts.at(pID).GetPosition(), parts.at(iPart).GetPosition(),
          parts.at(iPart).GetVelocity(), parts.at(iPart).GetAcceleration())};
      TVector3 pFieldv{pField.X(), pField.Y(), pField.Z()};
      particleField += pFieldv;
    }
  }
  return staticField + particleField;
}

TVector3 rad::ManyParticleSolver::acc(std::vector<Particle> &parts,
                                      size_t pID) {
  // Lorentz force
  TVector3 omega{get_omega(pID)};
  TVector3 acc{particles.at(pID).GetVelocity().Cross(omega)};

  // Add Larmor terms
  acc += radiation_acceleration(pID);

  // Acceleration from electric fields
  acc += particles.at(pID).GetCharge() * calc_e_field(parts, pID);

  return acc;
}

void rad::ManyParticleSolver::AdvanceParticles(double dt) {
  // Create a copy of the particles vector to calculate from
  std::vector<Particle> pInit{particles};

  for (size_t iPart{0}; iPart < particles.size(); iPart++) {
    Particle p{particles.at(iPart)};
    double gamma_n{1 / sqrt(1 - p.GetVelocity().Mag2() / pow(TMath::C(), 2))};
    TVector3 u_n{p.GetVelocity() * gamma_n};
    TVector3 x_n{p.GetPosition()};
    TVector3 v_n{p.GetVelocity()};

    // Half position step
    TVector3 x_nplushalf{x_n + v_n * (dt / 2.0)};

    // Do the first half of the Coloumb force
    TVector3 E_nplushalf{calc_e_field(pInit, iPart)};
    TVector3 E_tot_minus{E_nplushalf +
                         radiation_acceleration(x_nplushalf, u_n) *
                             (p.GetMass() / p.GetCharge())};
    TVector3 u_minus{u_n +
                     (dt * p.GetCharge() / (2 * p.GetMass())) * E_tot_minus};
    double gamma_minus{sqrt(1.0 + u_n.Mag2() / pow(TMath::C(), 2))};

    // Rotation step
    TVector3 B_nplushalf{calc_b_field(iPart)};
    double theta{p.GetCharge() * dt / (p.GetMass() * gamma_minus) *
                 B_nplushalf.Mag()};
    TVector3 u_minus_par{u_minus.Dot(B_nplushalf.Unit()) * B_nplushalf.Unit()};
    TVector3 u_plus{u_minus_par + (u_minus - u_minus_par) * cos(theta) +
                    (u_minus.Cross(B_nplushalf.Unit())) * sin(theta)};

    // Second half of the Coulomb force
    TVector3 E_tot_plus{E_nplushalf +
                        radiation_acceleration(x_nplushalf, u_plus) *
                            (p.GetMass() / p.GetCharge())};
    TVector3 u_nplus1{u_plus +
                      (dt * p.GetCharge() / (2 * p.GetMass())) * E_tot_plus};

    // Now update position
    double gamma_nplus1{sqrt(1 + pow(u_nplus1.Mag() / TMath::C(), 2))};
    TVector3 v_nplus1{u_nplus1 * (1 / gamma_nplus1)};
    TVector3 x_nplus1{x_nplushalf + v_nplus1 * (dt / 2.0)};

    // Update the particle position, velocity and acceleration
    particles.at(iPart).UpdateState(x_nplus1, v_nplus1, acc(pInit, iPart));
  }
}