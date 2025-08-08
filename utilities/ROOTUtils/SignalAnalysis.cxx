/// SignalAnalysis.cxx - ROOT-based signal analysis functions
#include "utilities/ROOTUtils/SignalAnalysis.h"

#include <cmath>
#include "utilities/BasicCore/Constants.h"

using namespace rad;

double rad::CalcAeHertzianDipole(const double wavelength,
                                 const ROOT::Math::XYZVector dipoleDir,
                                 const ROOT::Math::XYZPoint ePosition,
                                 const ROOT::Math::XYZPoint antennaPoint) {
  double Ae = 3 * pow(wavelength, 2) / (8 * PI);
  const double psi =
      std::acos(((ePosition - antennaPoint).Unit()).Dot(dipoleDir));
  Ae *= pow(std::sin(psi), 2);
  return Ae;
}

double rad::CalcAlHertzianDipole(const double wavelength,
                                 const ROOT::Math::XYZVector dipoleDir,
                                 const ROOT::Math::XYZPoint ePosition,
                                 const ROOT::Math::XYZPoint antennaPoint) {
  double Al = wavelength * std::sqrt(3 / (8 * PI));
  const double psi =
      std::acos(((ePosition - antennaPoint).Unit()).Dot(dipoleDir));
  Al *= std::sin(psi);
  return Al;
}

double rad::CalcRetardedTime(const ROOT::Math::XYZPoint fieldPoint,
                             const ROOT::Math::XYZPoint ePosition,
                             const double labTime) {
  double time =
      labTime - std::sqrt((ePosition - fieldPoint).Mag2()) / C;
  return time;
}

double rad::CalcTimeFromRetardedTime(ROOT::Math::XYZPoint fieldPoint,
                                     ROOT::Math::XYZPoint ePosition,
                                     double tRet) {
  double time =
      tRet + std::sqrt((ePosition - fieldPoint).Mag2()) / C;
  return time;
}

double rad::CalcTimeFromRetardedTime(TVector3 fieldPoint, TVector3 ePosition,
                                     double tRet) {
  double time = tRet + ((ePosition - fieldPoint).Mag() / C);
  return time;
}

double rad::GetGyroradius(TVector3 velocity, TVector3 bField,
                          double particleMass) {
  double gamma{1 /
               sqrt(1 - velocity.Dot(velocity) / (C * C))};
  TVector3 vPerp{velocity - (velocity.Dot(bField.Unit()) * bField)};
  double rg{gamma * particleMass * vPerp.Mag() / (QE * bField.Mag())};
  return rg;
}
