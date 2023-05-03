/*
  ElasticScatter.cxx
*/

#include "Scattering/ElasticScatter.h"

#include <cmath>

#include "TMath.h"

double rad::ElasticScatter::TotalRutherfordXSec() {
  double ke_keV{GetIncidentKE() / 1e3};
  double alpha{3.4e-3 / ke_keV};
  return 5.21e-25 * pow(1 / ke_keV, 2) * 4 * TMath::Pi() *
         pow((ke_keV + 511) / (ke_keV + 1022), 2) / (alpha * (1 + alpha));
}

double rad::ElasticScatter::GetTotalXSec() {
  double rutherfordxsec{TotalRutherfordXSec()};
  double gamma{0.734357 * (1 - exp(-2.45719 * sqrt(GetIncidentKE() / 1e3)))};
  return rutherfordxsec * gamma;
}