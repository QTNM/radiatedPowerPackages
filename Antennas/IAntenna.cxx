// IAntenna.cxx

#include "Antennas/IAntenna.h"

#include <cassert>
#include <iostream>

#include "BasicFunctions/BasicFunctions.h"
#include "TMath.h"
#include "TVector3.h"

double rad::IAntenna::GetCentralWavelength() {
  double freq = GetCentralFrequency();
  double lambda = TMath::C() / freq;
  return lambda;
}

TVector3 rad::IAntenna::GetThetaHat(const TVector3 electronPosition) {
  const TVector3 rHat = (antennaPosition - electronPosition).Unit();
  double x{antennaXAxis.Dot(rHat)};
  double y{antennaYAxis.Dot(rHat)};
  double z{antennaZAxis.Dot(rHat)};
  double theta{0};
  double phi{0};

  // Do theta first
  if (z > 0) {
    theta = atan(sqrt(x * x + y * y) / z);
  } else if (z < 0) {
    theta = atan(sqrt(x * x + y * y) / z) + TMath::Pi();
  } else {
    theta = TMath::Pi() / 2;
  }

  // Now do phi
  if (x > 0) {
    phi = atan(y / x);
  } else if (x < 0 && y >= 0) {
    phi = atan(y / x) + TMath::Pi();
  } else if (x < 0 && y < 0) {
    phi = atan(y / x) - TMath::Pi();
  } else if (x == 0 && y > 0) {
    phi = TMath::Pi() / 2;
  } else {
    phi = -TMath::Pi() / 2;
  }

  TVector3 vec(TMath::Cos(theta) * TMath::Cos(phi),
               TMath::Cos(theta) * TMath::Sin(phi), -1 * TMath::Sin(theta));
  TVector3 rotateVec{
      RotateToGlobalCoords(vec, antennaXAxis, antennaYAxis, antennaZAxis)};
  return rotateVec;
}

TVector3 rad::IAntenna::GetPhiHat(const TVector3 electronPosition) {
  const TVector3 rHat = (antennaPosition - electronPosition).Unit();
  double x{antennaXAxis.Dot(rHat)};
  double y{antennaYAxis.Dot(rHat)};
  double phi{0};
  if (x > 0) {
    phi = atan(y / x);
  } else if (x < 0 && y >= 0) {
    phi = atan(y / x) + TMath::Pi();
  } else if (x < 0 && y < 0) {
    phi = atan(y / x) - TMath::Pi();
  } else if (x == 0 && y > 0) {
    phi = TMath::Pi() / 2;
  } else {
    phi = -TMath::Pi() / 2;
  }

  TVector3 vec(-1 * TMath::Sin(phi), TMath::Cos(phi), 0);
  TVector3 rotateVec{
      RotateToGlobalCoords(vec, antennaXAxis, antennaYAxis, antennaZAxis)};
  return rotateVec;
}

double rad::IAntenna::GetTheta(const TVector3 electronPosition) {
  const TVector3 rHat = (antennaPosition - electronPosition).Unit();
  double theta = TMath::ACos(antennaZAxis.Dot(rHat));
  return theta;
}

double rad::IAntenna::GetPhi(const TVector3 electronPosition) {
  const TVector3 rHat = (antennaPosition - electronPosition).Unit();
  double phi = TMath::ATan2(antennaYAxis.Dot(rHat), antennaXAxis.Dot(rHat));
  if (phi < 0) phi += 2 * TMath::Pi();
  return phi;
}

void rad::IAntenna::SetBandwidth(const double lowerLimit,
                                 const double upperLimit) {
  assert(lowerLimit < upperLimit);
  lowerBandwidth = lowerLimit;
  upperBandwidth = upperLimit;
}

double rad::IAntenna::GetPatternIntegral() {
  // Now calculate the surface integral of the radiation pattern
  // Need this for proper calculation of the gain
  const int nPntsTheta{200};
  const int nPntsPhi{200};
  // Hypothetical width for integration
  const double binArea{TMath::Pi() * 2 * TMath::Pi() /
                       double(nPntsPhi * nPntsTheta)};
  const double binWidthTheta{TMath::Pi() / double(nPntsTheta)};
  const double binWidthPhi{2 * TMath::Pi() / double(nPntsPhi)};
  double PRad{0};
  for (int ith{0}; ith < nPntsTheta; ith++) {
    double theta{double(ith) * binWidthTheta + binWidthTheta / 2};
    for (int iph{0}; iph < nPntsPhi; iph++) {
      double phi{double(iph) * binWidthPhi + binWidthPhi / 2};
      double uSin{(GetETheta(theta, phi) * GetETheta(theta, phi) +
                   GetEPhi(theta, phi) * GetEPhi(theta, phi)) *
                  sin(theta)};
      PRad += uSin * binArea;
    }
  }
  return PRad;
}
