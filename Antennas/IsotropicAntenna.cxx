/// IsotropicAntenna.cxx

#include "Antennas/IsotropicAntenna.h"

#include <iostream>

#include "TVector3.h"

rad::IsotropicAntenna::IsotropicAntenna(TVector3 antPos, double AEff, double Z,
                                        double freq, double delay)
    : effectiveArea(AEff), imp(Z) {
  SetAntennaPosition(antPos);
  // Assume that antennas just point in world coordinate directions
  SetAntennaXAx(TVector3(1, 0, 0));
  SetAntennaYAx(TVector3(0, 1, 0));
  SetAntennaZAx(TVector3(0, 0, 1));
  SetTimeDelay(delay);
  SetBandwidth();
  SetCentralFreq(freq);

  PRad = GetPatternIntegral();
  std::cout << "PRad = " << PRad << std::endl;
}

TVector3 rad::IsotropicAntenna::GetETheta(const TVector3 ePos) {
  return GetThetaHat(ePos) * (1.0 / sqrt(2));
}

TVector3 rad::IsotropicAntenna::GetEPhi(const TVector3 ePos) {
  return GetPhiHat(ePos) * (1.0 / sqrt(2));
}

double rad::IsotropicAntenna::GetETheta(double theta, double phi) {
  return 1.0 / sqrt(2);
}

double rad::IsotropicAntenna::GetEPhi(double theta, double phi) {
  return 1.0 / sqrt(2);
}

double rad::IsotropicAntenna::GetHEff() {
  return sqrt(imp * effectiveArea / (120 * TMath::Pi()));
}

double rad::IsotropicAntenna::GetHEff(TVector3 ePos) { return GetHEff(); }

double rad::IsotropicAntenna::GetAEff(TVector3 electronPosition) {
  return effectiveArea;
}

double rad::IsotropicAntenna::GetAEffTheta(TVector3 ePos) {
  return effectiveArea / 2;
}

double rad::IsotropicAntenna::GetAEffPhi(TVector3 ePos) {
  return effectiveArea / 2;
}