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
    double rho, double phi, double z, Mode_t modeType, double A, unsigned int m,
    unsigned int n, unsigned int p, double t) {
  const double f{GetResonantModeF(modeType, m, n, p)};
  const double omega{f * TMath::TwoPi()};
  const double z0{d / 2};
  if (modeType == CircularCavity::kTE) {
    const double XmnPrime{GetBesselPrimeZero(m, n)};
    std::complex<double> Ez{0};
    double prefactor{A * omega / TMath::C()};
    prefactor *= pow(a / XmnPrime, 2);
    prefactor *= sin(double(p) * TMath::PiOver2() * (z / z0 + 1));

    std::complex<double> ERho1{
        -(double(m) / rho) * boost::math::cyl_bessel_j(m, XmnPrime * rho / a) *
        cos(omega * t - double(m) * phi)};
    ERho1 *= prefactor;
    std::complex<double> ERho2{
        (double(m) / rho) * boost::math::cyl_bessel_j(m, XmnPrime * rho / a) *
        cos(omega * t + double(m) * phi)};
    ERho2 *= prefactor;
    std::complex<double> EPhi1{
        (-XmnPrime / a) *
        boost::math::cyl_bessel_j_prime(m, XmnPrime * rho / a) *
        sin(omega * t - double(m) * phi)};
    EPhi1 *= prefactor;
    std::complex<double> EPhi2{
        (-XmnPrime / a) *
        boost::math::cyl_bessel_j_prime(m, XmnPrime * rho / a) *
        sin(omega * t + double(m) * phi)};
    EPhi2 *= prefactor;

    ComplexVector3 E1(ERho1 * cos(phi) - EPhi1 * sin(phi),
                      ERho1 * sin(phi) + EPhi1 * cos(phi), Ez);
    ComplexVector3 E2(ERho2 * cos(phi) - EPhi2 * sin(phi),
                      ERho2 * sin(phi) + EPhi2 * cos(phi), Ez);
    return /*E1+ */ E2;
  } else if (modeType == CircularCavity::kTM) {
    const double Xmn{boost::math::cyl_bessel_j_zero(double(m), n)};
    double prefactor{-double(p) * TMath::PiOver2() / z0};
    prefactor *=
        pow(a / Xmn, 2) * sin(double(p) * TMath::PiOver2() * (z / z0 + 1));
    std::complex<double> EPhi1{(double(m) / rho) *
                               boost::math::cyl_bessel_j(m, Xmn * rho / a) *
                               sin(omega * t - double(m) * phi)};
    EPhi1 *= prefactor * A;
    std::complex<double> EPhi2{(-double(m) / rho) *
                               boost::math::cyl_bessel_j(m, Xmn * rho / a) *
                               sin(omega * t + double(m) * phi)};
    EPhi2 *= prefactor * A;
    std::complex<double> ERho1{
        (Xmn / a) * boost::math::cyl_bessel_j_prime(m, Xmn * rho / a) *
        cos(omega * t - double(m) * phi)};
    ERho1 *= prefactor * A;
    std::complex<double> ERho2{
        (Xmn / a) * boost::math::cyl_bessel_j_prime(m, Xmn * rho / a) *
        cos(omega * t + double(m) * phi)};
    ERho2 *= prefactor * A;
    std::complex<double> Ez1{boost::math::cyl_bessel_j(m, Xmn * rho / a) *
                             cos(double(p) * TMath::PiOver2() * (z / z0 + 1)) *
                             cos(omega * t - double(m) * phi)};
    Ez1 *= A;
    std::complex<double> Ez2{boost::math::cyl_bessel_j(m, Xmn * rho / a) *
                             cos(double(p) * TMath::PiOver2() * (z / z0 + 1)) *
                             cos(omega * t + double(m) * phi)};
    Ez2 *= A;
    ComplexVector3 E1(ERho1 * cos(phi) - EPhi1 * sin(phi),
                      ERho1 * sin(phi) + EPhi1 * cos(phi), Ez1);
    ComplexVector3 E2(ERho2 * cos(phi) - EPhi2 * sin(phi),
                      ERho2 * sin(phi) + EPhi2 * cos(phi), Ez2);
    return E1;
  } else {
    std::cout << "TEM modes not supported for circular cavities!\n";
    return ComplexVector3(0, 0, 0);
  }
}

rad::ComplexVector3 rad::CircularCavity::GetModeEField(
    TVector3 pos, Mode_t modeType, double A, unsigned int n, unsigned int m,
    unsigned int l, double t) {
  double rho{pos.Perp()};
  double phi{atan2(pos.Y(), pos.X())};
  double z{pos.Z()};
  return GetModeEField(rho, phi, z, modeType, A, n, m, l, t);
}

rad::ComplexVector3 rad::CircularCavity::GetModeHField(
    double rho, double phi, double z, Mode_t modeType, double A, unsigned int n,
    unsigned int m, unsigned int l) {
  std::complex<double> j{0, 1};
  const double eta{376.730313668};
  const double f{GetResonantModeF(modeType, n, m, l)};
  const double k{TMath::TwoPi() * f / TMath::C()};

  // Check we are inside the cavity
  if (rho > a || z < 0 || z > d) {
    return ComplexVector3(0, 0, 0);
  } else {
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
}

rad::ComplexVector3 rad::CircularCavity::GetModeHField(TVector3 pos,
                                                       Mode_t modeType,
                                                       double A, unsigned int n,
                                                       unsigned int m,
                                                       unsigned int l) {
  double rho{pos.Perp()};
  double phi{atan2(pos.Y(), pos.X())};
  double z{pos.Z() + d / 2};
  return GetModeHField(rho, phi, z, modeType, A, n, m, l);
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