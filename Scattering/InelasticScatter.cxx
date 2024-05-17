/*
  InelasticScatter.cxx
*/

#include "Scattering/InelasticScatter.h"

#include <cmath>
#include <random>

#include "BasicFunctions/Constants.h"

double rad::InelasticScatter::GetTotalXSec() {
  double T{GetIncidentKE()};
  return S() * F(T) * g1(T);
}

double rad::InelasticScatter::G2(double W, double T) {
  double w{W / TRITIUM_I};
  double t{T / TRITIUM_I};
  return sqrt((w + 1) / t);
}

double rad::InelasticScatter::G3(double W, double T) {
  double w{W / TRITIUM_I};
  double t{T / TRITIUM_I};
  return BETA() * sqrt((1 - G2(W, T) * G2(W, T)) / w);
}

double rad::InelasticScatter::G4(double W, double T) {
  double w{W / TRITIUM_I};
  double t{T / TRITIUM_I};
  return GAMMA() * pow(1 - w / t, 3) / (t * (w + 1));
}

double rad::InelasticScatter::f_BE(double W, double T, double theta) {
  double w{W / TRITIUM_I};
  double t{T / TRITIUM_I};
  return 1 / (1 + pow((cos(theta) - G2(W, T)) / G3(W, T), 2));
}

double rad::InelasticScatter::f_b(double theta) {
  return 1.0 / (1 + pow((cos(theta) + 1) / G5(), 2));
}

double rad::InelasticScatter::g_BE(double W, double T) {
  double term1{atan((1 - G2(W, T)) / G3(W, T))};
  double term2{atan((1 + G2(W, T)) / G3(W, T))};
  return 2 * M_PI * G3(W, T) * (term1 + term2);
}

double rad::InelasticScatter::F(double T) {
  double t{T / TRITIUM_I};
  return (A1() * log(t) + A2() + A3() / t) / t;
}

double rad::InelasticScatter::f_1(double W, double T) {
  double w{W / TRITIUM_I};
  double t{T / TRITIUM_I};
  double term1{1 / pow(w + 1, n())};
  double term2{1 / pow(t - w, n())};
  double term3{-1 / pow((w + 1) * (t - w), n() / 2)};
  return term1 + term2 + term3;
}

double rad::InelasticScatter::G1(double W, double T) {
  return (S() * F(T) * f_1(W, T) / TRITIUM_I) / (g_BE(W, T) + G4(W, T) * G_B());
}

double rad::InelasticScatter::G4fb(double W, double T, double theta) {
  return G4(W, T) / (1 + pow((cos(theta) + 1) / G5(), 2));
}

double rad::InelasticScatter::g1(double T) {
  double t{T / TRITIUM_I};
  return (1 - pow(t, 1 - n())) / (n() - 1) -
         pow(2 / (t + 1), n() / 2) * (1 - pow(t, 1 - n() / 2)) / (n() - 2);
}

double rad::InelasticScatter::GetDoubleDiffXSec(double W, double theta) {
  double T{GetIncidentKE()};
  return G1(W, T) * (f_BE(W, T, theta) + G4fb(W, T, theta));
}

double rad::InelasticScatter::GetSDCS_W(double W) {
  double T{GetIncidentKE()};
  return G1(W, T) * (g_BE(W, T) + G4(W, T) * G_B());
}

double rad::InelasticScatter::GetRandomW() {
  const unsigned int nBins{400};
  // We want the differences between the points to be logarithmically spaced
  const double diffMin{1e-5};
  const double diffMax{GetIncidentKE() / 2};
  const double WDiff{log10(diffMax / diffMin) / double(nBins - 1)};

  // Build the distribution of differential cross-section
  std::vector<double> WVec{};
  std::vector<double> xsecVec{};
  for (int iW{nBins - 1}; iW >= 0; iW--) {
    double W{GetIncidentKE() - TRITIUM_I -
             diffMin * pow(10, double(iW) * WDiff)};
    WVec.push_back(W);
    xsecVec.push_back(GetSDCS_W(W));
  }

  // Sample from the distribution
  std::piecewise_linear_distribution<double> diffXSecDist(
      WVec.begin(), WVec.end(), xsecVec.begin());
  std::random_device rd;
  std::mt19937 mt(rd());
  return diffXSecDist(mt);
}

double rad::InelasticScatter::GetRandomTheta(double W) {
  const unsigned int nBins{1000};
  const double thetaMin{0};
  const double thetaMax{M_PI};
  const double thetaBinWidth{(thetaMax - thetaMin) / double(nBins)};

  std::vector<double> thetaVec{};
  std::vector<double> xsecVec{};
  for (unsigned int iT{0}; iT < nBins; iT++) {
    double theta{thetaMin + thetaBinWidth / 2 + double(iT) * thetaBinWidth};
    thetaVec.push_back(theta);
    xsecVec.push_back(GetDoubleDiffXSec(W, theta));
  }

  // Sample from the distribution
  std::piecewise_linear_distribution<double> xsecDist(
      thetaVec.begin(), thetaVec.end(), xsecVec.begin());
  std::random_device rd;
  std::mt19937 mt(rd());
  return xsecDist(mt);
}

double rad::InelasticScatter::GetPrimaryScatteredE(double W, double theta) {
  double E1{GetIncidentKE()};
  double E2Prime{W};
  return E1 - E2Prime - RYDBERG_EV;
}

double rad::InelasticScatter::GetPrimaryScatteredAngle(double W, double theta) {
  double E1{GetIncidentKE()};
  double E1Prime{GetPrimaryScatteredE(W, theta)};
  double E2Prime{W};
  double cTheta1{(E1 + E1Prime - E2Prime) / (2 * sqrt(E1 * E1Prime))};
  return acos(cTheta1);
}