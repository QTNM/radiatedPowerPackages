// IAntenna.cxx

#include "physics/Antennas/IAntenna.h"

#include <cassert>
#include <cmath>
#include <iostream>

#include "utilities/BasicFunctions/BasicFunctions.h"
#include "utilities/BasicCore/Constants.h"
#include "TVector3.h"

using namespace rad;

double rad::IAntenna::GetCentralWavelength() {
  double freq = GetCentralFrequency();
  double lambda = C / freq;
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
    theta = atan(sqrt(x * x + y * y) / z) + PI;
  } else {
    theta = PI / 2;
  }

  // Now do phi
  if (x > 0) {
    phi = atan(y / x);
  } else if (x < 0 && y >= 0) {
    phi = atan(y / x) + PI;
  } else if (x < 0 && y < 0) {
    phi = atan(y / x) - PI;
  } else if (x == 0 && y > 0) {
    phi = PI / 2;
  } else {
    phi = -PI / 2;
  }

  TVector3 vec(std::cos(theta) * std::cos(phi),
               std::cos(theta) * std::sin(phi), -1 * std::sin(theta));
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
    phi = atan(y / x) + PI;
  } else if (x < 0 && y < 0) {
    phi = atan(y / x) - PI;
  } else if (x == 0 && y > 0) {
    phi = PI / 2;
  } else {
    phi = -PI / 2;
  }

  TVector3 vec(-1 * std::sin(phi), std::cos(phi), 0);
  TVector3 rotateVec{
      RotateToGlobalCoords(vec, antennaXAxis, antennaYAxis, antennaZAxis)};
  return rotateVec;
}

double rad::IAntenna::GetTheta(const TVector3 electronPosition) {
  const TVector3 rHat = (antennaPosition - electronPosition).Unit();
  double theta = std::acos(antennaZAxis.Dot(rHat));
  return theta;
}

double rad::IAntenna::GetPhi(const TVector3 electronPosition) {
  const TVector3 rHat = (antennaPosition - electronPosition).Unit();
  double phi = std::atan2(antennaYAxis.Dot(rHat), antennaXAxis.Dot(rHat));
  if (phi < 0) phi += 2 * PI;
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
  const double binArea{PI * 2 * PI /
                       double(nPntsPhi * nPntsTheta)};
  const double binWidthTheta{PI / double(nPntsTheta)};
  const double binWidthPhi{2 * PI / double(nPntsPhi)};
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
