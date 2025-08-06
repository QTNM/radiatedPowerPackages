/// TritiumSpectrum.cxx

#include "BasicFunctions/TritiumSpectrum.h"

#include "BasicCore/Constants.h"
#include "BasicCore/NuFitValues.h"
#include "BasicFunctions/BasicFunctions.h"
#include "TMath.h"

double rad::TritiumDecayRate(double electronT, double m1, double m2, double m3,
                             double endpointE) {
  double beta = sqrt(1.0 - pow(ME_EV / (electronT + ME_EV), 2));
  double p = sqrt(electronT * electronT + 2 * ME_EV * electronT);
  double eta = ALPHA * 2 / beta;
  double fermiFunction =
      2 * TMath::Pi() * eta / (1 - exp(-2 * TMath::Pi() * eta));
  double nuE = endpointE - electronT;
  double rate = G_F * G_F * pow(0.97425, 2) * fermiFunction *
                (1 + 3.0 * pow(-1.2646, 2)) * p * (electronT + ME_EV) /
                (2 * pow(TMath::Pi(), 3));

  // PMNS matrix elements
  double Ue1Sq = pow(cos(kNuFitTh12NH) * cos(kNuFitTh13NH), 2);
  double Ue2Sq = pow(sin(kNuFitTh12NH) * cos(kNuFitTh13NH), 2);
  double Ue3Sq = pow(sin(kNuFitTh13NH), 2);

  double neutrinoPhaseSpc1 =
      (m1 <= nuE) ? Ue1Sq * nuE * sqrt(nuE * nuE - m1 * m1) : 0.0;
  double neutrinoPhaseSpc2 =
      (m2 <= nuE) ? Ue2Sq * nuE * sqrt(nuE * nuE - m2 * m2) : 0.0;
  double neutrinoPhaseSpc3 =
      (m3 <= nuE) ? Ue3Sq * nuE * sqrt(nuE * nuE - m3 * m3) : 0.0;
  double neutrinoPhaseSpc =
      neutrinoPhaseSpc1 + neutrinoPhaseSpc2 + neutrinoPhaseSpc3;
  rate *= neutrinoPhaseSpc;
  rate *= 1.0 / 6.58e-16;  // Account for natural units

  return rate;
}
