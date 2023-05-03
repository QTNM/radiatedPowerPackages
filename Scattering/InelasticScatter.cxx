/*
  InelasticScatter.cxx
*/

#include "Scattering/InelasticScatter.h"

#include <cmath>

#include "BasicFunctions/Constants.h"

double rad::InelasticScatter::GetTotalXSec() {
  double U{GetIncidentKE() / RYDBERG_EV};
  return 1.3 * TMath::Pi() * A0 * A0 * (log(U) + 4) / U;
}