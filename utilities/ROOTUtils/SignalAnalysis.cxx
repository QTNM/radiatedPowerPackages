/// SignalAnalysis.cxx - ROOT-based signal analysis functions
#include "utilities/ROOTUtils/SignalAnalysis.h"

#include "utilities/BasicCore/Constants.h"
#include "TMath.h"

double rad::CalcAeHertzianDipole(const double wavelength,
                                 const ROOT::Math::XYZVector dipoleDir,
                                 const ROOT::Math::XYZPoint ePosition,
                                 const ROOT::Math::XYZPoint antennaPoint) {
  double Ae = 3 * pow(wavelength, 2) / (8 * TMath::Pi());
  const double psi =
      TMath::ACos(((ePosition - antennaPoint).Unit()).Dot(dipoleDir));
  Ae *= pow(TMath::Sin(psi), 2);
  return Ae;
}

double rad::CalcAlHertzianDipole(const double wavelength,
                                 const ROOT::Math::XYZVector dipoleDir,
                                 const ROOT::Math::XYZPoint ePosition,
                                 const ROOT::Math::XYZPoint antennaPoint) {
  double Al = wavelength * TMath::Sqrt(3 / (8 * TMath::Pi()));
  const double psi =
      TMath::ACos(((ePosition - antennaPoint).Unit()).Dot(dipoleDir));
  Al *= TMath::Sin(psi);
  return Al;
}

double rad::CalcRetardedTime(const ROOT::Math::XYZPoint fieldPoint,
                             const ROOT::Math::XYZPoint ePosition,
                             const double labTime) {
  double time =
      labTime - TMath::Sqrt((ePosition - fieldPoint).Mag2()) / TMath::C();
  return time;
}

double rad::CalcTimeFromRetardedTime(ROOT::Math::XYZPoint fieldPoint,
                                     ROOT::Math::XYZPoint ePosition,
                                     double tRet) {
  double time =
      tRet + TMath::Sqrt((ePosition - fieldPoint).Mag2()) / TMath::C();
  return time;
}

double rad::CalcTimeFromRetardedTime(TVector3 fieldPoint, TVector3 ePosition,
                                     double tRet) {
  double time = tRet + ((ePosition - fieldPoint).Mag() / TMath::C());
  return time;
}

double rad::GetSpeedFromKE(double T, double particleMass) {
  double gamma = T * TMath::Qe() / (ME * TMath::C() * TMath::C()) + 1;
  double betaSq = 1 - 1 / pow(gamma, 2);
  double speed = sqrt(betaSq) * TMath::C();
  return speed;
}

double rad::GetGyroradius(TVector3 velocity, TVector3 bField,
                          double particleMass) {
  double gamma{1 /
               sqrt(1 - velocity.Dot(velocity) / (TMath::C() * TMath::C()))};
  TVector3 vPerp{velocity - (velocity.Dot(bField.Unit()) * bField)};
  double rg{gamma * particleMass * vPerp.Mag() / (TMath::Qe() * bField.Mag())};
  return rg;
}

TVector3 rad::calculate_omega(const TVector3 BField, const double charge,
                              const double energy, const double mass) {
  double gamma_m0 = mass + energy * TMath::Qe() / pow(TMath::C(), 2);
  return (charge * BField * (1.0 / gamma_m0));
}

double rad::CalcCyclotronFreq(const double KE, const double B) {
  double freq =
      TMath::Qe() * B / (ME + (KE * TMath::Qe() / pow(TMath::C(), 2)));
  return freq / (2 * TMath::Pi());
}