// PatchAntenna.cxx

#include "Antennas/PatchAntenna.h"

#include <cassert>
#include <iostream>

#include "TMath.h"
#include "TVector3.h"
#include "boost/math/special_functions/sinc.hpp"

rad::PatchAntenna::PatchAntenna(TVector3 antPos, TVector3 antXAx,
                                TVector3 antYAx, double width, double height,
                                double f0, double permittivity, double delay) {
  antYAx = antYAx.Unit();
  antXAx = antXAx.Unit();

  // Make sure axes are perpendicular
  assert(antXAx.Dot(antYAx) == 0);

  SetAntennaPosition(antPos);
  SetAntennaXAx(antXAx);
  SetAntennaYAx(antYAx);
  SetAntennaZAx(antXAx.Cross(antYAx));

  H = height;
  W = width;

  relativePerm = permittivity;
  SetCentralFreq(f0);
  SetTimeDelay(delay);

  // First calculate the effective perimittivity
  double epEff{((permittivity + 1.0) / 2.0) +
               ((permittivity - 1.0) / 2.0) * pow(1.0 + 12.0 * H / W, -0.5)};
  // Then calculate the effective length
  LEff = TMath::C() / (2 * f0 * sqrt(epEff));
  // From this, we can calculate the actual length of the patch
  double deltaL{0.412 * H * (epEff + 0.3) * (W / H + 0.264) /
                ((epEff - 0.258) * (W / H + 0.8))};
  L = LEff - 2 * deltaL;

  SetBandwidth();

  PRad = GetPatternIntegral();
}

TVector3 rad::PatchAntenna::GetETheta(const TVector3 electronPosition) {
  TVector3 thetaHat = GetThetaHat(electronPosition);
  double thetaAng = GetTheta(electronPosition);
  double phiAng = GetPhi(electronPosition);
  return thetaHat * GetETheta(thetaAng, phiAng);
}

TVector3 rad::PatchAntenna::GetEPhi(const TVector3 electronPosition) {
  TVector3 phiHat = GetPhiHat(electronPosition);
  double thetaAng = GetTheta(electronPosition);
  double phiAng = GetPhi(electronPosition);
  return phiHat * GetEPhi(thetaAng, phiAng);
}

double rad::PatchAntenna::GetETheta(double theta, double phi) {
  if (theta > TMath::Pi() / 2) {
    return 0;
  } else {
    double k{2 * TMath::Pi() * GetCentralFreq() / TMath::C()};
    double F1{boost::math::sinc_pi(k * H * sin(theta) * cos(theta) / 2.0) *
              boost::math::sinc_pi(k * W * sin(theta) * sin(phi) / 2.0)};
    double F2{2.0 * cos(k * L * sin(theta) * cos(phi) / 2.0)};
    double ETheta{cos(phi) * F1 * F2};
    return ETheta;
  }
}

double rad::PatchAntenna::GetEPhi(double theta, double phi) {
  if (theta > TMath::Pi() / 2) {
    return 0;
  } else {
    double k{2 * TMath::Pi() * GetCentralFreq() / TMath::C()};
    double F1{boost::math::sinc_pi(k * H * sin(theta) * cos(theta) / 2.0) *
              boost::math::sinc_pi(k * W * sin(theta) * sin(phi) / 2.0)};
    double F2{2.0 * cos(k * L * sin(theta) * cos(phi) / 2.0)};
    double EPhi{cos(theta) * sin(phi) * F1 * F2};
    return EPhi;
  }
}

double rad::PatchAntenna::GetHEff(TVector3 ePos) {
  double theta{GetTheta(ePos)};
  double phi{GetPhi(ePos)};
  double gain{4 * TMath::Pi() *
              (GetETheta(theta, phi) * GetETheta(theta, phi) +
               GetEPhi(theta, phi) * GetEPhi(theta, phi)) /
              PRad};
  double imp{50};  // Assume we can match to 50 Ohm load
  double lambda{TMath::C() / GetCentralFreq()};
  return sqrt(imp * lambda * lambda * gain / (480 * TMath::Pi() * TMath::Pi()));
}

double rad::PatchAntenna::GetAEff(TVector3 ePos) {
  double theta{GetTheta(ePos)};
  double phi{GetPhi(ePos)};
  double gain{4 * TMath::Pi() *
              (GetETheta(theta, phi) * GetETheta(theta, phi) +
               GetEPhi(theta, phi) * GetEPhi(theta, phi)) /
              PRad};
  return pow(TMath::C() / GetCentralFreq(), 2) * gain / (4 * TMath::Pi());
}

double rad::PatchAntenna::GetImpedance() {
  double Z{90 * (pow(relativePerm, 2) / (relativePerm - 1)) * pow(L / W, 2)};
  return Z;
}

double rad::PatchAntenna::GetBandwidth() {
  double lambda{TMath::C() / GetCentralFreq()};
  double bw{GetCentralFreq() * 3.77 * W * H *
            (relativePerm - 1 / pow(relativePerm, 2)) / (L * lambda)};
  return bw;
}

double rad::PatchAntenna::GetAEffTheta(TVector3 ePos) {
  double theta{GetTheta(ePos)};
  double phi{GetPhi(ePos)};
  double gain{4 * TMath::Pi() * GetETheta(theta, phi) * GetETheta(theta, phi) /
              PRad};
  return pow(TMath::C() / GetCentralFreq(), 2) * gain / (4 * TMath::Pi());
}

double rad::PatchAntenna::GetAEffPhi(TVector3 ePos) {
  double theta{GetTheta(ePos)};
  double phi{GetPhi(ePos)};
  double gain{4 * TMath::Pi() * GetEPhi(theta, phi) * GetEPhi(theta, phi) /
              PRad};
  return pow(TMath::C() / GetCentralFreq(), 2) * gain / (4 * TMath::Pi());
}