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

double rad::InelasticScatter::G2(double omega, double t) {
  return sqrt((omega + 1) / t);
}

double rad::InelasticScatter::G3(double omega, double t) {
  double g2{G2(omega, t)};
  return BETA() * sqrt((1 - g2 * g2) / omega);
}

double rad::InelasticScatter::G4(double omega, double t) {
  return GAMMA() * pow(1 - omega / t, 3) / (t * (omega + 1));
}

double rad::InelasticScatter::f_BE(double omega, double t, double theta) {
  double g2{G2(omega, t)};
  double g3{G3(omega, t)};
  return 1.0 / (1 + pow((cos(theta) - g2) / g3, 2));
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

double rad::InelasticScatter::GetDoubleDiffXSec(double W, double theta) {
  double omega{W / RYDBERG_EV};
  double t{GetIncidentKE() / RYDBERG_EV};
  double g1{G1(omega, t)};
  double g4{G4(omega, t)};
  return g1 * (f_BE(omega, t, theta) + g4 * f_b(theta));
}