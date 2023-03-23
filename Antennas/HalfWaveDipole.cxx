// HalfWaveDipole.cxx

#include "Antennas/HalfWaveDipole.h"

#include <cassert>

#include "TVector3.h"

rad::HalfWaveDipole::HalfWaveDipole(TVector3 antPos, TVector3 antXAx,
                                    TVector3 antZAx, double freq,
                                    double delay) {
  // Make sure all axes are unit vectors initially
  antXAx = antXAx.Unit();
  antZAx = antZAx.Unit();

  // Make sure that axes are perpendicular to one another
  assert(antXAx.Dot(antZAx) == 0);

  SetAntennaPosition(antPos);
  SetAntennaXAx(antXAx);
  SetAntennaYAx(antZAx.Cross(antXAx));
  SetAntennaZAx(antZAx);
  SetCentralFreq(freq);
  SetTimeDelay(delay);

  SetBandwidth();

  PRad = GetPatternIntegral();
}

// Calculate the radiation pattern in the theta hat direction
TVector3 rad::HalfWaveDipole::GetETheta(const TVector3 electronPosition) {
  TVector3 thetaHat = GetThetaHat(electronPosition);
  double thetaAng = GetTheta(electronPosition);
  thetaHat *=
      TMath::Cos(TMath::Pi() * TMath::Cos(thetaAng) / 2) / TMath::Sin(thetaAng);
  return thetaHat;
}

double rad::HalfWaveDipole::GetETheta(double theta, double phi) {
  return cos(TMath::Pi() * cos(theta) / 2) / sin(theta);
}

double rad::HalfWaveDipole::GetEPhi(double theta, double phi) { return 0; }

TVector3 rad::HalfWaveDipole::GetEPhi(const TVector3 electronPosition) {
  // No radiation in the phi direction for a hertzian dipole
  return TVector3(0, 0, 0);
}

double rad::HalfWaveDipole::GetHEff() {
  double heff = GetCentralWavelength() / TMath::Pi();
  return heff;
}

double rad::HalfWaveDipole::GetHEff(TVector3 ePos) {
  double theta{GetTheta(ePos)};
  double phi{GetPhi(ePos)};
  double gain{4 * TMath::Pi() *
              (GetETheta(theta, phi) * GetETheta(theta, phi) +
               GetEPhi(theta, phi) * GetEPhi(theta, phi)) /
              PRad};
  double imp{73};  // Characteristic impedance of half-wave dipole
  double lambda{TMath::C() / GetCentralFreq()};
  return sqrt(imp * lambda * lambda * gain / (480 * TMath::Pi() * TMath::Pi()));
}

double rad::HalfWaveDipole::GetAEff(TVector3 ePos) {
  // Gain of a half-wave dipole is 1.65 at theta = pi / 2
  double theta{GetTheta(ePos)};
  double phi{GetPhi(ePos)};
  double gain{4 * TMath::Pi() *
              (GetETheta(theta, phi) * GetETheta(theta, phi) +
               GetEPhi(theta, phi) * GetEPhi(theta, phi)) /
              PRad};
  return pow(TMath::C() / GetCentralFreq(), 2) * gain / (4 * TMath::Pi());
}

double rad::HalfWaveDipole::GetAEffTheta(TVector3 ePos) {
  double theta{GetTheta(ePos)};
  double phi{GetPhi(ePos)};
  double gain{4 * TMath::Pi() * GetETheta(theta, phi) * GetETheta(theta, phi) /
              PRad};
  return pow(TMath::C() / GetCentralFreq(), 2) * gain / (4 * TMath::Pi());
}

double rad::HalfWaveDipole::GetAEffPhi(TVector3 ePos) { return 0; }