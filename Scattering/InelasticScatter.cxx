/*
  InelasticScatter.cxx
*/

#include "Scattering/InelasticScatter.h"

#include <cmath>
#include <random>

#include "BasicFunctions/Constants.h"

inline double rad::InelasticScatter::BETA() { return 0.60; }

inline double rad::InelasticScatter::GAMMA() { return 10.0; }

inline double rad::InelasticScatter::G_B() { return 2.9; }

inline double rad::InelasticScatter::n() { return 2.5; }

inline double rad::InelasticScatter::A1() { return 0.74; }

inline double rad::InelasticScatter::A2() { return 0.87; }

inline double rad::InelasticScatter::A3() { return -0.6; }

inline double rad::InelasticScatter::S() {
  return 4 * TMath::Pi() * pow(A0, 2);
}

double rad::InelasticScatter::GetTotalXSec() {
  double U{GetIncidentKE() / RYDBERG_EV};
  return 1.3 * TMath::Pi() * A0 * A0 * (log(U) + 4) / U;
}

double rad::InelasticScatter::G2(double omega, double t) {
  return sqrt((omega + 1) / t);
}

double rad::InelasticScatter::G3(double omega, double t) {
  return BETA() * sqrt((1 - pow(G2(omega, t), 2)) / omega);
}

double rad::InelasticScatter::G4(double omega, double t) {
  return GAMMA() * pow(1 - omega / t, 3) / (t * (omega + 1));
}

inline double rad::InelasticScatter::G5() { return 0.33; }

double rad::InelasticScatter::f_BE(double omega, double t, double theta) {
  return 1.0 / (1 + pow((cos(theta) - G2(omega, t)) / G3(omega, t), 2));
}

double rad::InelasticScatter::f_b(double theta) {
  return 1.0 / (1 + pow((cos(theta) + 1) / G5(), 2));
}

double rad::InelasticScatter::g_BE(double omega, double t) {
  double g2{G2(omega, t)};
  double g3{G3(omega, t)};
  double term1{atan2((1 - g2), g3)};
  double term2{atan2((1 + g2), g3)};
  return 2 * TMath::Pi() * g3 * (term1 + term2);
}

double rad::InelasticScatter::F(double t) {
  return (A1() * log(t) + A2() + A3() / t) / t;
}

double rad::InelasticScatter::f_1(double omega, double t) {
  double term1{1 / pow(omega + 1, n())};
  double term2{1 / pow(t - omega, n())};
  double term3{1 / pow((omega + 1) * (t - omega), n() / 2)};
  return term1 + term2 + term3;
}

double rad::InelasticScatter::G1(double omega, double t) {
  double numerator{S() * F(t) * f_1(omega, t) / RYDBERG_EV};
  double denom{g_BE(omega, t) + G4(omega, t) * G_B()};
  return numerator / denom;
}

double rad::InelasticScatter::GetSingleDiffXSec_W(double W) {
  const double omega{W / RYDBERG_EV};
  const double t{GetIncidentKE() / RYDBERG_EV};
  return G1(omega, t) * (g_BE(omega, t) + G4(omega, t) * G_B());
}

double rad::InelasticScatter::GetDoubleDiffXSec(double W, double theta) {
  double omega{W / RYDBERG_EV};
  double t{GetIncidentKE() / RYDBERG_EV};
  double g1{G1(omega, t)};
  double g4{G4(omega, t)};
  return g1 * (f_BE(omega, t, theta) + g4 * f_b(theta));
}

double rad::InelasticScatter::GetRandomW() {
  const unsigned int nBins{200};
  const double omegaMin{0};
  const double omegaMax{8};
  const double wMin{omegaMin * RYDBERG_EV};
  const double wMax{omegaMax * RYDBERG_EV};
  const double wBinWidth{(wMax - wMin) / double(nBins)};
  const double thetaMin{0};
  const double thetaMax{TMath::PiOver2()};
  const double thetaBinWidth{(thetaMax - thetaMin) / double(nBins)};

  // Build the distribution of differential cross-section
  std::vector<double> wVec{};
  std::vector<double> xsecVec{};
  for (unsigned int iW{0}; iW < nBins; iW++) {
    double diffXSec{0};
    double w{wMin + wBinWidth / 2 + double(iW) * wBinWidth};
    wVec.push_back(w);
    for (unsigned int iT{0}; iT < nBins; iT++) {
      double theta{thetaMin + thetaBinWidth / 2 + double(iT) * thetaBinWidth};
      diffXSec += GetDoubleDiffXSec(w, theta) * wBinWidth;
    }
    xsecVec.push_back(diffXSec);
  }

  // Sample from the distribution
  std::piecewise_linear_distribution<double> diffXSecDist(
      wVec.begin(), wVec.end(), xsecVec.begin());
  std::random_device rd;
  std::mt19937 mt(rd());
  return diffXSecDist(mt);
}

double rad::InelasticScatter::GetRandomTheta(double W) {
  const unsigned int nBins{200};
  const double thetaMin{0};
  const double thetaMax{TMath::PiOver2()};
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