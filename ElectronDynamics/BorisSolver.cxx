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
                              const double mass_v, const double tau_v,
                              CircularCavity* cavity)
    : charge(charge_v), mass(mass_v), field(field_v), tau(tau_v), cav(cavity) {
  // If a cavity is present then calculate the Purcell factor as we go
  if (cav && tau != 0) {
    calcPurcellFactor = true;

    // Firstly need to calculate field normalisations for the cavity modes
    // Normalise so that field = 1 at antinode
    const unsigned int nScanPnts{50};  // Scan 50 points in each direction
    const double A{1};
    for (unsigned int iRho{0}; iRho < nScanPnts; iRho++) {
      double rho{double(iRho) * cav->GetRadius() / double(nScanPnts - 1)};
      for (unsigned int iPhi{0}; iPhi < nScanPnts; iPhi++) {
        double phi{double(iPhi) * TMath::TwoPi() / double(nScanPnts - 1)};
        for (unsigned int iZ{0}; iZ < nScanPnts; iZ++) {
          double z{-cav->GetLength() +
                   cav->GetLength() * double(iZ) / double(nScanPnts - 1)};
          TVector3 pos(rho * cos(phi), rho * sin(phi), z);
          // For now just look at the TE111 mode
          TVector3 fPlus{
              cav->GetModalEField(pos, CircularCavity::kTE, A, 1, 1, 1, true)
                  .Real()};
          TVector3 fMinus{
              cav->GetModalEField(pos, CircularCavity::kTE, A, 1, 1, 1, false)
                  .Real()};

          if (fPlus.Mag() > maxFieldPlus) maxFieldPlus = fPlus.Mag();
          if (fMinus.Mag() > maxFieldMinus) maxFieldMinus = fMinus.Mag();
          if ((fPlus + fMinus).Mag() > maxField)
            maxField = (fPlus + fMinus).Mag();
        }
      }
    }
    std::cout << "Field maximum +ve, -ve = " << maxFieldPlus << ", "
              << maxFieldMinus << std::endl;
    std::cout << "Combined field maximum = " << maxField << std::endl;

    // Now calculate the effective volume of the cavity
    // Integrate field over cavity volume
    const double dRho{cav->GetRadius() / double(nScanPnts)};
    const double dPhi{TMath::TwoPi() / double(nScanPnts)};
    const double dZ{cav->GetLength() / double(nScanPnts)};
    for (unsigned int iRho{0}; iRho < nScanPnts; iRho++) {
      double rho{dRho / 2 + double(iRho) * dRho};
      double dV{rho * dRho * dPhi * dZ};
      for (unsigned int iPhi{0}; iPhi < nScanPnts; iPhi++) {
        double phi{dPhi / 2 + double(iPhi) * dPhi};
        for (unsigned int iZ{0}; iZ < nScanPnts; iZ++) {
          double z{dZ / 2 + double(iZ) * dZ};
          TVector3 pos(rho * cos(phi), rho * sin(phi), z);

          TVector3 fPlus{cav->GetModalEField(pos, CircularCavity::kTE,
                                             1 / maxFieldPlus, 1, 1, 1, true)
                             .Real()};
          vEffPlus += dV * fPlus.Mag();
          TVector3 fMinus{cav->GetModalEField(pos, CircularCavity::kTE,
                                              1 / maxFieldMinus, 1, 1, 1, false)
                              .Real()};
          vEffMinus += dV * fMinus.Mag();
          vEff += dV * (fPlus + fMinus).Mag();
        }
      }
    }
    std::cout << "Effective volume +ve, -ve = " << vEffPlus * 1e9 << " mm^3, "
              << vEffMinus * 1e9 << " mm^3\n";
    std::cout << "Effective volume = " << vEff * 1e9 << " mm^3\n";

    // Get resonant frequency
    const double cavQ{200};
    const double fTE111{cav->GetResonantModeF(CircularCavity::kTE, 1, 1, 1)};
    const double lTE111{TMath::C() / fTE111};
    FpMax = 3 * cavQ * pow(lTE111, 3) / (pow(TMath::TwoPi(), 2) * vEff);
    std::cout << "Maximum Purcell factor = " << FpMax << std::endl;
  }
}

TVector3 rad::BorisSolver::get_omega(const TVector3 pos) {
  TVector3 BField = field->evaluate_field_at_point(pos);
  return calculate_omega(BField, charge, 0.0, mass);
}

TVector3 rad::BorisSolver::radiation_acceleration(const TVector3 pos,
                                                  const TVector3 vel) {
  double fieldFactor{1};
  double detuningFactor{1};
  if (cav && tau != 0) {
    // If a cavity is present then do some calculation of the relative mode
    // field strength
    ComplexVector3 fPlus{cav->GetModalEField(pos, CircularCavity::kTE,
                                             1 / maxField, 1, 1, 1, true)};
    ComplexVector3 fMinus{cav->GetModalEField(pos, CircularCavity::kTE,
                                              1 / maxField, 1, 1, 1, false)};
    fieldFactor = (fPlus + fMinus).Real().Mag2();

    // Now calculate the detuning factor
    // Start by getting the mode resonant frequency
    double fRes{cav->GetResonantModeF(CircularCavity::kTE, 1, 1, 1)};
    // Now get the instantaneous magnetic field at this point
    double gamma{1 / sqrt(1 - pow(vel.Mag() / TMath::C(), 2))};
    double ke{(gamma - 1) * ME * TMath::C() * TMath::C() / TMath::Qe()};
    double f{CalcCyclotronFreq(ke, calc_b_field(pos).Mag())};
    double deltaFRes{fRes / 200};

    // Calculate detuning based upon how far we are from the resonance
    detuningFactor =
        deltaFRes * deltaFRes / (4 * pow(f - fRes, 2) + deltaFRes * deltaFRes);
  }

  double Fp{FpMax * fieldFactor * detuningFactor};

  TVector3 omega{get_omega(pos)};
  double denom{1 + tau * tau * omega.Dot(omega)};
  double accX{0};
  double accY{0};
  double accZ{0};

  accX -= tau * (omega.Z() * omega.Z() + omega.Y() * omega.Y()) * vel.X();
  accX += tau * omega.X() * (omega.Z() * vel.Z() + omega.Y() * vel.Y());

  accY -= tau * (omega.Z() * omega.Z() + omega.X() * omega.X()) * vel.Y();
  accY += tau * omega.Y() * (omega.Z() * vel.Z() + omega.X() * vel.X());

  accZ -= tau * (omega.X() * omega.X() + omega.Y() * omega.Y()) * vel.Z();
  accZ += tau * omega.Z() * (omega.X() * vel.X() + omega.Y() * vel.Y());

  TVector3 acc(accX / denom, accY / denom, accZ / denom);
  return acc * sqrt(Fp);
}

TVector3 rad::BorisSolver::acc(const TVector3 pos, const TVector3 vel) {
  TVector3 omega{get_omega(pos)};

  // Lorentz force
  TVector3 acc{vel.Cross(omega)};

  // Add Larmor terms
  acc += radiation_acceleration(pos, vel);

  // Acceleration from electric field
  acc += charge * calc_e_field(pos);

  return acc;
}

std::tuple<TVector3, TVector3> rad::BorisSolver::advance_step(
    const double time_step, const TVector3 x0, const TVector3 v0) {
  double gamma_n{1.0 / sqrt(1 - v0.Dot(v0) / pow(TMath::C(), 2))};
  TVector3 u_n{v0 * gamma_n};
  TVector3 x_n{x0};
  TVector3 v_n{v0};

  // Half position step
  TVector3 x_nplushalf{x_n + v_n * (time_step / 2.0)};

  // Do the first half of the Coloumb force
  TVector3 E_nplushalf{calc_e_field(x_nplushalf)};
  TVector3 E_tot_minus{E_nplushalf + radiation_acceleration(x_nplushalf, u_n) *
                                         (mass / charge)};
  TVector3 u_minus{u_n + (time_step * charge / (2 * mass)) * E_tot_minus};
  double gamma_minus{sqrt(1.0 + u_n.Dot(u_n) / pow(TMath::C(), 2))};

  // Rotation step
  TVector3 B_nplushalf{calc_b_field(x_nplushalf)};
  double theta{charge * time_step / (mass * gamma_minus) * B_nplushalf.Mag()};
  TVector3 u_minus_par{u_minus.Dot(B_nplushalf.Unit()) * B_nplushalf.Unit()};
  TVector3 u_plus{u_minus_par + (u_minus - u_minus_par) * cos(theta) +
                  (u_minus.Cross(B_nplushalf.Unit())) * sin(theta)};

  // Second half of the Coulomb force
  TVector3 E_tot_plus{E_nplushalf +
                      radiation_acceleration(x_nplushalf, u_plus) *
                          (mass / charge)};
  TVector3 u_nplus1{u_plus + (time_step * charge / (2 * mass)) * E_tot_plus};

  // Now update position
  double gamma_nplus1{sqrt(1 + pow(u_nplus1.Mag() / TMath::C(), 2))};
  TVector3 v_nplus1{u_nplus1 * (1 / gamma_nplus1)};
  TVector3 x_nplus1{x_nplushalf + v_nplus1 * (time_step / 2.0)};

  return std::make_tuple(x_nplus1, v_nplus1);
}

TVector3 rad::BorisSolver::calc_b_field(TVector3 pos) {
  return (field->evaluate_field_at_point(pos));
}

TVector3 rad::BorisSolver::calc_e_field(TVector3 pos) {
  return (field->evaluate_e_field_at_point(pos));
}