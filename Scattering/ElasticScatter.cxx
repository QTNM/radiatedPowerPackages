/*
  ElasticScatter.cxx
*/

#include "Scattering/ElasticScatter.h"

#include <boost/math/special_functions/legendre.hpp>
#include <boost/units/systems/si/codata/universal_constants.hpp>
#include <cmath>
#include <iostream>
#include <random>

#include "BasicFunctions/BasicFunctions.h"
#include "TMath.h"

rad::ElasticScatter::ElasticScatter(double incidentKE, unsigned int aNum,
                                    unsigned int aMass)
    : BaseScatter(incidentKE), Z(aNum), A(aMass) {
  // Check if we have data for this atomic number
  if (Z > 5) {
    std::cout << "No data for this atomic number, returning zero!" << std::endl;
  }
}

double rad::ElasticScatter::TotalRutherfordXSec() {
  double ke_keV{GetIncidentKE() / 1e3};
  double alpha{3.4e-3 * pow(Z, 2 / 3) / ke_keV};
  return 5.21e-25 * pow(Z / ke_keV, 2) * 4 * TMath::Pi() *
         pow((ke_keV + ME_EV / 1e3) / (ke_keV + 2 * ME_EV / 1e3), 2) /
         (alpha * (1 + alpha));
}

double rad::ElasticScatter::RutherfordDCS(double theta) {
  const double chargeFac{
      Z * pow(TMath::Qe(), 2) /
      pow(16 * M_PI * EPSILON0 * GetIncidentKE() * TMath::Qe(), 2)};
  if (theta == 0) {
    return 0;
  } else {
    return chargeFac / pow(sin(theta / 2), 4);
  }
}

double rad::ElasticScatter::Calcq(double theta) {
  const double gamma{1 + GetIncidentKE() / ME_EV};
  const double beta{sqrt(1 - 1 / pow(gamma, 2))};
  const double p{gamma * beta * ME * TMath::C()};
  return 2 * p * sin(theta / 2);
}

double rad::ElasticScatter::Calct() {
  const double gamma{1 + GetIncidentKE() / ME_EV};
  const double beta{sqrt(1 - 1 / pow(gamma, 2))};
  return 1 + 1.053 * pow(ALPHA * Z, 1.971) +
         (0.00569 + 0.995 * pow(ALPHA * Z, 1.778) * pow(1 / beta - 1, 0.750));
}

double rad::ElasticScatter::CalcFe(double q) {
  double Fe{0};
  for (size_t i{0}; i < 4; i++) {
    Fe += Ai[Z - 1][i] * pow(alphai[Z - 1][i], 2) /
          (pow(alphai[Z - 1][i], 2) + pow(q / 1.992851915e-24, 2));
  }
  return Fe;
}

double rad::ElasticScatter::ScreeningFactor(double theta) {
  return pow(1 - CalcFe(Calcq(theta) / Calct()), 2);
}

double rad::ElasticScatter::CalcF(double R, double q) {
  return 3 / pow(q * R / TMath::Hbar(), 3) *
         (sin(q * R / TMath::Hbar()) -
          q * R / TMath::Hbar() * cos(q * R / TMath::Hbar()));
}

inline double rad::ElasticScatter::R0() { return 1.2e-15 * pow(A, 1 / 3); }

double rad::ElasticScatter::NuclearFormFactor(double q) {
  return CalcF(R0(), 2) * CalcF(2e-15, q);
}

double rad::ElasticScatter::GetTotalXSec() {
  double rutherfordxsec{TotalRutherfordXSec()};
  double lambda{0};
  double beta{0};
  if (Z == 1) {
    lambda = 0.734357;
    beta = 2.45719;
  } else if (Z == 2) {
    lambda = 1.25;
    beta = 12;
  } else if (Z == 3) {
    lambda = 1.262598;
    beta = 10.17333;
  } else if (Z == 4) {
    lambda = 1.492947;
    beta = 4.38619;
  } else if (Z == 5) {
    lambda = 1.434886;
    beta = 4.016957;
  }
  double gamma{lambda * (1 - exp(-beta * sqrt(GetIncidentKE() / 1e3)))};
  return rutherfordxsec * gamma;
}

double rad::ElasticScatter::GetDiffXSec(std::vector<double> A, double B,
                                        std::vector<double> C,
                                        double cosTheta) {
  double sum1{0};
  double sum2{0};
  for (int m{1}; m <= 4; m++) {
    sum1 += A.at(m - 1) * pow(1 - cosTheta + 2 * B, -m);
  }

  for (unsigned int n{0}; n <= 6; n++) {
    sum2 += C.at(n) * boost::math::legendre_p(n, cosTheta);
  }

  return sum1 + sum2;
}

double rad::ElasticScatter::GetDiffXSec(double theta) {
  if (Z > 5) {
    return 0;
  } else {
    const double q{Calcq(theta)};
    return RutherfordDCS(theta) * ScreeningFactor(theta) * NuclearFormFactor(q);
  }
}

double rad::ElasticScatter::GetRandomScatteringAngle() {
  const unsigned int nPnts{1000};
  const double thetaMin{1e-5};
  const double thetaMax{179 * M_PI /
                        180};  // Strange behaviour close to 180 degrees
  // Calculate the step size in log space
  const double stepSize{log10(thetaMax / thetaMin) / (nPnts - 1)};

  std::vector<double> thetaVec{};
  std::vector<double> dcsVec{};
  for (unsigned int i{0}; i < nPnts; i++) {
    double theta{thetaMin * pow(10, i * stepSize)};
    thetaVec.push_back(theta);
    double dcs{GetDiffXSec(theta)};
    dcsVec.push_back(dcs);
  }

  // Now construct the distribution and sample from it
  std::piecewise_linear_distribution<double> dcsDist(
      thetaVec.begin(), thetaVec.end(), dcsVec.begin());
  std::random_device rd;
  std::mt19937 mt(rd());
  return dcsDist(mt);
}

double rad::ElasticScatter::GetEnergyAfterScatter(double theta) {
  // Calculate incident momentum
  const double Ei_Joules{GetIncidentKE() * TMath::Qe() +
                         ME * pow(TMath::C(), 2)};
  const double pi{sqrt(pow(Ei_Joules, 2) - pow(ME * pow(TMath::C(), 2), 2)) /
                  TMath::C()};
  // Now calculate the final momentum
  const double pf{pi - Calcq(theta)};
  const double Ef{
      sqrt(pow(pf * TMath::C(), 2) + pow(ME * pow(TMath::C(), 2), 2))};
  return (Ef - ME * pow(TMath::C(), 2)) / TMath::Qe();
}