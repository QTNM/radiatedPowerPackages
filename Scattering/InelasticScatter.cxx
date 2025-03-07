/*
  InelasticScatter.cxx
*/

#include "Scattering/InelasticScatter.h"

#include <cmath>
#include <iostream>
#include <random>

#include "BasicFunctions/Constants.h"

rad::InelasticScatter::InelasticScatter(double T, Species species)
    : BaseScatter(T), theSpecies(species) {}

double rad::InelasticScatter::I() {
  if (theSpecies == H2) {
    return 15.43;
  } else if (theSpecies == H) {
    return RYDBERG_EV;
  } else if (theSpecies == He) {
    return 24.59;
  } else {
    std::cerr << "Species not supported! Returning 0.\n";
    return 0;
  }
}

double rad::InelasticScatter::U() {
  if (theSpecies == H2) {
    return 25.68;
  } else if (theSpecies == H) {
    return RYDBERG_EV;
  } else {
    return 39.51;
  }
}

double rad::InelasticScatter::S() {
  return 4 * M_PI * N() * pow(A0 * RYDBERG_EV / I(), 2);
}

double rad::InelasticScatter::N() {
  if (theSpecies == H2) {
    return 2;
  } else if (theSpecies == H) {
    return 1;
  } else {
    return 2;
  }
}

double rad::InelasticScatter::Ni() {
  if (theSpecies == H2) {
    return 1.173;
  } else if (theSpecies == H) {
    return 0.4343;
  } else {
    return 1.605;
  }
}

double rad::InelasticScatter::D() {
  const double t{GetIncidentKE() / I()};
  double b{0}, c{0}, d{0}, e{0}, f{0};
  if (theSpecies == H2) {
    c = 1.1262;
    d = 6.3982;
    e = -7.8055;
    f = 2.144;
  } else if (theSpecies == H) {
    b = -2.2473e-2;
    c = 1.1775;
    d = -4.6264e-1;
    e = 8.9064e-2;
  } else {
    c = 1.2178e1;
    d = -2.9585e1;
    e = 3.1251e1;
    f = -1.2175e1;
  }
  double tTerm{(t + 1) / 2};
  double bTerm{(b / 2) * (1 - pow(tTerm, -2))};
  double cTerm{(c / 3) * (1 - pow(tTerm, -3))};
  double dTerm{(d / 4) * (1 - pow(tTerm, -4))};
  double eTerm{(e / 5) * (1 - pow(tTerm, -5))};
  double fTerm{(f / 6) * (1 - pow(tTerm, -6))};
  return (bTerm + cTerm + dTerm + eTerm + fTerm) / N();
}

double rad::InelasticScatter::GetTotalXSec() {
  const double t{GetIncidentKE() / I()};
  const double u{U() / I()};
  return S() / (t + u + 1) *
         (D() * log(t) + (2 - Ni() / N()) * ((t - 1) / t - log(t) / (t + 1)));
}

inline double rad::InelasticScatter::BETA() { return 0.60; }

inline double rad::InelasticScatter::GAMMA() { return 10.0; }

inline double rad::InelasticScatter::G_B() { return 2.9; }

inline double rad::InelasticScatter::n() { return 2.4; }

double rad::InelasticScatter::A1() {
  if (theSpecies == He) {
    return 0.85;
  } else {
    return 0.74;
  }
}

double rad::InelasticScatter::A2() {
  if (theSpecies == He) {
    return 0.36;
  } else {
    return 0.87;
  }
}

double rad::InelasticScatter::A3() {
  if (theSpecies == He) {
    return -0.10;
  } else {
    return -0.60;
  }
}

double rad::InelasticScatter::G2(double W, double T) {
  double w{W / I()};
  double t{T / I()};
  return sqrt((w + 1) / t);
}

double rad::InelasticScatter::G3(double W, double T) {
  double w{W / I()};
  double t{T / I()};
  return BETA() * sqrt((1 - G2(W, T) * G2(W, T)) / w);
}

double rad::InelasticScatter::G4(double W, double T) {
  double w{W / I()};
  double t{T / I()};
  return GAMMA() * pow(1 - w / t, 3) / (t * (w + 1));
}

inline double rad::InelasticScatter::G5() { return 0.33; }

double rad::InelasticScatter::f_BE(double W, double T, double theta) {
  double w{W / I()};
  double t{T / I()};
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
  double t{T / I()};
  return (A1() * log(t) + A2() + A3() / t) / t;
}

double rad::InelasticScatter::f_1(double W, double T) {
  double w{W / I()};
  double t{T / I()};
  double term1{1 / pow(w + 1, n())};
  double term2{1 / pow(t - w, n())};
  double term3{-1 / pow((w + 1) * (t - w), n() / 2)};
  return term1 + term2 + term3;
}

double rad::InelasticScatter::G1(double W, double T) {
  return (S() * F(T) * f_1(W, T) / I()) / (g_BE(W, T) + G4(W, T) * G_B());
}

double rad::InelasticScatter::G4fb(double W, double T, double theta) {
  return G4(W, T) / (1 + pow((cos(theta) + 1) / G5(), 2));
}

double rad::InelasticScatter::g1(double T) {
  double t{T / I()};
  return (1 - pow(t, 1 - n())) / (n() - 1) -
         pow(2 / (t + 1), n() / 2) * (1 - pow(t, 1 - n() / 2)) / (n() - 2);
}

double rad::InelasticScatter::GetDoubleDiffXSec(double W, double theta) {
  double T{GetIncidentKE()};
  return G1(W, T) * (f_BE(W, T, theta) + G4fb(W, T, theta));
}

double rad::InelasticScatter::GetSDCS_W(double W) {
  const double t{GetIncidentKE() / I()};
  const double u{U() / I()};
  const double w{W / I()};
  double b{0}, c{0}, d{0}, e{0}, f{0};
  if (theSpecies == H2) {
    c = 1.1262;
    d = 6.3982;
    e = -7.8055;
    f = 2.144;
  } else if (theSpecies == H) {
    b = -2.2473e-2;
    c = 1.1775;
    d = -4.6264e-1;
    e = 8.9064e-2;
  } else {
    c = 1.2178e1;
    d = -2.9585e1;
    e = 3.1251e1;
    f = -1.2175e1;
  }
  const double diffOsc{b * pow(w + 1, -2) + c * pow(w + 1, -3) +
                       d * pow(w + 1, -4) + e * pow(w + 1, -5) +
                       f * pow(w + 1, -6)};
  const double preFactor{S() / (I() * (t + u + 1))};
  const double term1{(Ni() / N() - 2) * (1 / (w + 1) + 1 / (t - w)) / (t + 1)};
  const double term2{(2 - Ni() / N()) * (pow(w + 1, -2) + pow(t - w, -2))};
  const double term3{log(t) * diffOsc / (N() * (w + 1))};
  return preFactor * (term1 + term2 + term3);
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
    double W{GetIncidentKE() - I() - diffMin * pow(10, double(iW) * WDiff)};
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