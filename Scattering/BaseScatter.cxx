/*
  BaseScatter.cxx

  Class implementation
*/

#include "Scattering/BaseScatter.h"

#include "BasicFunctions/BasicFunctions.h"
#include "BasicFunctions/Constants.h"

double rad::BaseScatter::GetMeanFreePath(double N) {
  double xsec{GetTotalXSec()};
  return 1 / (N * xsec);
}

double rad::BaseScatter::GetMeanFreeTime(double N) {
  double lambda{GetMeanFreePath(N)};
  double v{GetSpeedFromKE(ke, ME)};
  return lambda / v;
}