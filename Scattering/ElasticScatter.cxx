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
  double alpha{3.4e-3 / ke_keV};
  return 5.21e-25 * pow(1 / ke_keV, 2) * 4 * TMath::Pi() *
         pow((ke_keV + 511) / (ke_keV + 1022), 2) / (alpha * (1 + alpha));
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
  double gamma{0.734357 * (1 - exp(-2.45719 * sqrt(GetIncidentKE() / 1e3)))};
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

double rad::ElasticScatter::GetDiffXSec(double cosTheta) {
  double angstrom{1e-10};  // metres
  double T{GetIncidentKE()};

  // Figure out which elements to use for the interpolation
  unsigned int iLo{0};
  if (T >= 8e3 && T < 16e3) {
    iLo = 1;
  } else if (T >= 16e3) {
    iLo = 2;
  }

  std::vector<double> eTmp{eVec.at(iLo), eVec.at(iLo + 1), eVec.at(iLo + 2),
                           eVec.at(iLo + 3)};
  std::vector<double> bTmp{BVec.at(iLo), BVec.at(iLo + 1), BVec.at(iLo + 2),
                           BVec.at(iLo + 3)};
  std::vector<double> a1Tmp{A1Vec.at(iLo), A1Vec.at(iLo + 1), A1Vec.at(iLo + 2),
                            A1Vec.at(iLo + 3)};
  std::vector<double> a2Tmp{A2Vec.at(iLo), A2Vec.at(iLo + 1), A2Vec.at(iLo + 2),
                            A2Vec.at(iLo + 3)};
  std::vector<double> a3Tmp{A3Vec.at(iLo), A3Vec.at(iLo + 1), A3Vec.at(iLo + 2),
                            A3Vec.at(iLo + 3)};
  std::vector<double> a4Tmp{A4Vec.at(iLo), A4Vec.at(iLo + 1), A4Vec.at(iLo + 2),
                            A4Vec.at(iLo + 3)};
  std::vector<double> c0Tmp{C0Vec.at(iLo), C0Vec.at(iLo + 1), C0Vec.at(iLo + 2),
                            C0Vec.at(iLo + 3)};
  std::vector<double> c1Tmp{C1Vec.at(iLo), C1Vec.at(iLo + 1), C1Vec.at(iLo + 2),
                            C1Vec.at(iLo + 3)};
  std::vector<double> c2Tmp{C2Vec.at(iLo), C2Vec.at(iLo + 1), C2Vec.at(iLo + 2),
                            C2Vec.at(iLo + 3)};
  std::vector<double> c3Tmp{C3Vec.at(iLo), C3Vec.at(iLo + 1), C3Vec.at(iLo + 2),
                            C3Vec.at(iLo + 3)};
  std::vector<double> c4Tmp{C4Vec.at(iLo), C4Vec.at(iLo + 1), C4Vec.at(iLo + 2),
                            C4Vec.at(iLo + 3)};
  std::vector<double> c5Tmp{C5Vec.at(iLo), C5Vec.at(iLo + 1), C5Vec.at(iLo + 2),
                            C5Vec.at(iLo + 3)};
  std::vector<double> c6Tmp{C6Vec.at(iLo), C6Vec.at(iLo + 1), C6Vec.at(iLo + 2),
                            C6Vec.at(iLo + 3)};

  std::vector<double> aVec{};
  aVec.push_back(CubicInterpolation(eTmp, a1Tmp, T));
  aVec.push_back(CubicInterpolation(eTmp, a2Tmp, T));
  aVec.push_back(CubicInterpolation(eTmp, a3Tmp, T));
  aVec.push_back(CubicInterpolation(eTmp, a4Tmp, T));
  double b{CubicInterpolation(eTmp, bTmp, T)};
  std::vector<double> cVec{};
  cVec.push_back(CubicInterpolation(eTmp, c0Tmp, T));
  cVec.push_back(CubicInterpolation(eTmp, c1Tmp, T));
  cVec.push_back(CubicInterpolation(eTmp, c2Tmp, T));
  cVec.push_back(CubicInterpolation(eTmp, c3Tmp, T));
  cVec.push_back(CubicInterpolation(eTmp, c4Tmp, T));
  cVec.push_back(CubicInterpolation(eTmp, c5Tmp, T));
  cVec.push_back(CubicInterpolation(eTmp, c6Tmp, T));

  return GetDiffXSec(aVec, b, cVec, cosTheta) * angstrom * angstrom;
}

double rad::ElasticScatter::GetRandomScatteringAngle() {
  const unsigned int nPnts{1000};
  const double cThetaMin{-1};
  const double cThetaMax{0.9995};  // Strange behaviour close to 1

  std::vector<double> cThetaVec{};
  std::vector<double> xsecVec{};
  for (unsigned int i{0}; i < nPnts; i++) {
    double cTheta{cThetaMin +
                  (cThetaMax - cThetaMin) * double(i) / double(nPnts - 1)};
    cThetaVec.push_back(cTheta);
    double xsec{GetDiffXSec(cTheta)};
    xsecVec.push_back(xsec);
  }

  // Now construct the distribution and sample from it
  std::piecewise_linear_distribution<double> xsecDist(
      cThetaVec.begin(), cThetaVec.end(), xsecVec.begin());
  std::random_device rd;
  std::mt19937 mt(rd());
  double cThetaSampled{xsecDist(mt)};
  return acos(cThetaSampled);
}