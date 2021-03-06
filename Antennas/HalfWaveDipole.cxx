// HalfWaveDipole.cxx

#include <cassert>

#include "Antennas/HalfWaveDipole.h"

#include "TVector3.h"

rad::HalfWaveDipole::HalfWaveDipole(TVector3 antPos, TVector3 antXAx, TVector3 antZAx,
				    double freq, double delay) {
  // Make sure all axes are unit vectors initially
  antXAx = antXAx.Unit();
  antZAx = antZAx.Unit();

  // Make sure that axes are perpendicular to one another
  assert(antXAx.Dot(antZAx) == 0);
  
  antennaPosition = antPos;
  antennaXAxis = antXAx;
  antennaYAxis = antZAx.Cross(antXAx);
  antennaZAxis = antZAx;
  centralFreq = freq;
  timeDelay = delay;
  
  SetBandwidth();
}

// Calculate the radiation pattern in the theta hat direction
TVector3 rad::HalfWaveDipole::GetETheta(const TVector3 electronPosition) {
  TVector3 thetaHat = GetThetaHat(electronPosition);
  double thetaAng = GetTheta(electronPosition);
  thetaHat *= TMath::Cos(TMath::Pi() * TMath::Cos(thetaAng) / 2) / TMath::Sin(thetaAng);
  return thetaHat;
}

TVector3 rad::HalfWaveDipole::GetEPhi(const TVector3 electronPosition) {
  // No radiation in the phi direction for a hertzian dipole
  return TVector3(0, 0, 0);
}

double rad::HalfWaveDipole::GetHEff() {
  double heff = GetCentralWavelength() / TMath::Pi();
  return heff;
}
