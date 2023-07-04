#include "Waveguides/CircularCavity.h"

#include <boost/math/special_functions/bessel.hpp>
#include <boost/math/special_functions/bessel_prime.hpp>
#include <complex>
#include <iostream>

#include "BasicFunctions/BasicFunctions.h"
#include "TMath.h"

double rad::CircularCavity::GetResonantModeF(Mode_t modeType, unsigned int n,
                                             unsigned int m, unsigned int l) {
  if (m < 1) {
    std::cout << "m cannot be 0! Returning -1.\n";
    return -1;
  } else {
    if (modeType == CircularCavity::kTE) {
      const double pnmPrime{GetBesselPrimeZero(n, m)};
      double fSq{pow(TMath::C() / TMath::TwoPi(), 2) *
                 (pow(pnmPrime / a, 2) + pow(double(l) * TMath::Pi() / d, 2))};
      return sqrt(fSq);
    } else if (modeType == CircularCavity::kTM) {
      double pnm{boost::math::cyl_bessel_j_zero(double(n), m)};
      double fSq{pow(TMath::C() / TMath::TwoPi(), 2) *
                 (pow(pnm / a, 2) + pow(double(l) * TMath::Pi() / d, 2))};
      return sqrt(fSq);
    } else {
      std::cout << "TEM modes not supported for circular cavities!\n";
      return -1;
    }
  }
}

rad::ComplexVector3 rad::CircularCavity::GetModeEField(
    double rho, double phi, double z, Mode_t modeType, double A, unsigned int n,
    unsigned int m, unsigned int l) {
  std::complex<double> j{0, 1};
  const double eta{376.730313668};
  const double f{GetResonantModeF(modeType, n, m, l)};
  const double k{TMath::TwoPi() * f / TMath::C()};

  if (modeType == CircularCavity::kTE) {
    const double pnmPrime{GetBesselPrimeZero(n, m)};
    std::complex<double> Ez{0};
    std::complex<double> ERho{j * k * eta * A * a * a * double(n) /
                              (pnmPrime * pnmPrime * rho)};
    ERho *= sin(double(n) * phi) * sin(double(l) * TMath::Pi() * z / d);
    ERho *= boost::math::cyl_bessel_j(n, pnmPrime * rho / a);
    std::complex<double> EPhi{j * k * eta * a * A / pnmPrime};
    EPhi *= cos(double(n) * phi) * sin(double(l) * TMath::Pi() * z / d);
    EPhi *= boost::math::cyl_bessel_j_prime(n, pnmPrime * rho / a);
    ComplexVector3 E(ERho * cos(phi) - EPhi * sin(phi),
                     ERho * sin(phi) + EPhi * cos(phi), Ez);
    return E;
  } else {
    return ComplexVector3(0, 0, 0);
  }
}

rad::ComplexVector3 rad::CircularCavity::GetModeHField(
    double rho, double phi, double z, Mode_t modeType, double A, unsigned int n,
    unsigned int m, unsigned int l) {
  std::complex<double> j{0, 1};
  const double eta{376.730313668};
  const double f{GetResonantModeF(modeType, n, m, l)};
  const double k{TMath::TwoPi() * f / TMath::C()};

  if (modeType == CircularCavity::kTE) {
    const double pnmPrime{GetBesselPrimeZero(n, m)};
    std::complex<double> beta{sqrt(k * k - pow(pnmPrime / a, 2))};
    std::complex<double> Hz{A *
                            boost::math::cyl_bessel_j(n, pnmPrime * rho / a)};
    Hz *= cos(double(n) * phi) * sin(double(l) * TMath::Pi() * z / d);
    std::complex<double> HRho{beta * a * A / pnmPrime};
    HRho *= boost::math::cyl_bessel_j_prime(double(n), pnmPrime * rho / a);
    HRho *= cos(double(n) * phi) * cos(double(l) * TMath::Pi() * z / d);
    std::complex<double> HPhi{-beta * a * a * double(n) * A /
                              (pnmPrime * pnmPrime * rho)};
    HPhi *= boost::math::cyl_bessel_j(double(n), pnmPrime * rho / a);
    HPhi *= sin(double(n) * phi) * cos(double(l) * TMath::Pi() * z / d);
    ComplexVector3 H(HRho * cos(phi) - HPhi * sin(phi),
                     HRho * sin(phi) + HPhi * cos(phi), Hz);
    return H;
  } else {
    return ComplexVector3(0, 0, 0);
  }
}

double rad::CircularCavity::GetCutoffFrequency(Mode_t modeType, int n, int m) {
  if (modeType == kTE) {
    double pnmPrime{GetBesselPrimeZero(n, m)};
    double k_c{pnmPrime / a};
    double f_c{k_c * TMath::C() / (2 * TMath::Pi())};
    return f_c;
  } else if (modeType == kTM) {
    double pnm{boost::math::cyl_bessel_j_zero(double(n), m)};
    double k_c{pnm / a};
    double f_c{k_c * TMath::C() / (2 * TMath::Pi())};
    return f_c;
  } else {
    std::cout << "Unsupported mode type! Returning -1!\n";
    return -1;
  }
}