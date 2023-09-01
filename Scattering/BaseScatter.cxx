/*
  BaseScatter.cxx

  Class implementation
*/

#include "Scattering/BaseScatter.h"

#include <random>

#include "BasicFunctions/BasicFunctions.h"
#include "BasicFunctions/Constants.h"
#include "TMath.h"

double rad::BaseScatter::GetMeanFreePath(double N) {
  double xsec{GetTotalXSec()};
  return 1 / (N * xsec);
}

double rad::BaseScatter::GetMeanFreeTime(double N) {
  double lambda{GetMeanFreePath(N)};
  double v{GetSpeedFromKE(ke, ME)};
  return lambda / v;
}

TVector3 rad::BaseScatter::GetScatteredVector(TVector3 vel, double outKE,
                                              double theta) {
  // Set up RNG
  std::random_device rd;
  std::mt19937 gen(rd());

  // Calculate outgoing speed from kinetic energy
  const double speed{TMath::C() * sqrt((pow(ME_EV + outKE, 2) - ME_EV * ME_EV) /
                                       pow(ME_EV + outKE, 2))};

  // Distribute azimuthal angle uniformly
  std::uniform_real_distribution<double> uni(0, 1);
  double phi{uni(gen) * TMath::TwoPi()};

  // Original direction in global coords
  TVector3 originalDir{vel.Unit()};
  // Original direction (velocity along z axis)
  TVector3 oldDir(0, 0, 1);
  TVector3 newDir(sin(theta) * cos(phi), sin(theta) * sin(phi), cos(theta));
  // Define other axes
  TVector3 ax2(-originalDir.Y() / originalDir.X(), 1, 0);
  ax2 = ax2.Unit();
  TVector3 ax3{(ax2.Cross(originalDir)).Unit()};

  // Do the actual rotation
  TVector3 newDirection{RotateToCoords(newDir, ax2, ax3, originalDir)};
  return speed * newDirection;
}