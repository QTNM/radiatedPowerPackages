/*
  ElasticScatter.cxx
*/

#include "Scattering/ElasticScatter.h"

#include <boost/math/special_functions/legendre.hpp>
#include <cmath>
#include <random>

#include "BasicFunctions/BasicFunctions.h"
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