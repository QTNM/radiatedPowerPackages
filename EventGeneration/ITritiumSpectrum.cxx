/// ITritiumSpectrum.cxx

#include "EventGeneration/ITritiumSpectrum.h"

#include <cmath>

double rad::ITritiumSpectrum::CalculateUe1Sq()
{
  double ue1sq = pow( cos(th12)*cos(th13), 2 );
  return ue1sq;
}

double rad::ITritiumSpectrum::CalculateUe2Sq()
{
  double ue2sq = pow( sin(th12)*cos(th13), 2 );
  return ue2sq;
}

double rad::ITritiumSpectrum::CalculateUe3Sq()
{
  double ue3sq = pow( sin(th13), 2 );
  return ue3sq;
}
