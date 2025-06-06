#include "Waveguides/CircularCavity.h"

#include <boost/math/special_functions/bessel.hpp>
#include <boost/math/special_functions/bessel_prime.hpp>
#include <complex>
#include <iostream>

#include "BasicFunctions/BasicFunctions.h"
#include "TMath.h"

rad::CircularCavity::CircularCavity(double radius, double length,
                                    TVector3 probe)
    : a(radius), d(length) {
  SetProbePosition(probe);
}

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
    unsigned int n, unsigned int p, bool state, double t) {
  const double f{GetResonantModeF(modeType, m, n, p)};
  const double omega{f * TMath::TwoPi()};
  const double z0{d / 2};

  // If we're outside the cavity then the field should be 0
  if (rho > a || abs(z) > z0) return ComplexVector3(0, 0, 0);

  std::complex<double> EPhi{0};
  std::complex<double> ERho{0};
  std::complex<double> EZ{0};
  if (modeType == CircularCavity::kTE) {
    const double XmnPrime{GetBesselPrimeZero(m, n)};
    double prefactor{A * omega / TMath::C()};
    prefactor *= pow(a / XmnPrime, 2);
    prefactor *= sin(double(p) * TMath::PiOver2() * (z / z0 + 1));

    if (state) {
      ERho = -(double(m) / rho) *
             boost::math::cyl_bessel_j(m, XmnPrime * rho / a) *
             cos(omega * t - double(m) * phi);
      ERho *= prefactor;
      EPhi = (-XmnPrime / a) *
             boost::math::cyl_bessel_j_prime(m, XmnPrime * rho / a) *
             sin(omega * t - double(m) * phi);
      EPhi *= prefactor;
    } else {
      ERho = (double(m) / rho) *
             boost::math::cyl_bessel_j(m, XmnPrime * rho / a) *
             cos(omega * t + double(m) * phi);
      ERho *= prefactor;
      EPhi = (-XmnPrime / a) *
             boost::math::cyl_bessel_j_prime(m, XmnPrime * rho / a) *
             sin(omega * t + double(m) * phi);
      EPhi *= prefactor;
    }
  } else if (modeType == CircularCavity::kTM) {
    const double Xmn{boost::math::cyl_bessel_j_zero(double(m), n)};
    double prefactor{-double(p) * TMath::PiOver2() / z0};
    prefactor *=
        pow(a / Xmn, 2) * sin(double(p) * TMath::PiOver2() * (z / z0 + 1));

    if (state) {
      EPhi = (double(m) / rho) * boost::math::cyl_bessel_j(m, Xmn * rho / a) *
             sin(omega * t - double(m) * phi);
      EPhi *= prefactor * A;
      ERho = (Xmn / a) * boost::math::cyl_bessel_j_prime(m, Xmn * rho / a) *
             cos(omega * t - double(m) * phi);
      ERho *= prefactor * A;
      EZ = boost::math::cyl_bessel_j(m, Xmn * rho / a) *
           cos(double(p) * TMath::PiOver2() * (z / z0 + 1)) *
           cos(omega * t - double(m) * phi);
      EZ *= A;
    } else {
      EPhi = (-double(m) / rho) * boost::math::cyl_bessel_j(m, Xmn * rho / a) *
             sin(omega * t + double(m) * phi);
      EPhi *= prefactor * A;

      ERho = (Xmn / a) * boost::math::cyl_bessel_j_prime(m, Xmn * rho / a) *
             cos(omega * t + double(m) * phi);
      ERho *= prefactor * A;

      EZ = boost::math::cyl_bessel_j(m, Xmn * rho / a) *
           cos(double(p) * TMath::PiOver2() * (z / z0 + 1)) *
           cos(omega * t + double(m) * phi);
      EZ *= A;
    }
  } else {
    std::cout << "TEM modes not supported for circular cavities!\n";
  }

  // Generate final field vector
  ComplexVector3 E(ERho * cos(phi) - EPhi * sin(phi),
                   ERho * sin(phi) + EPhi * cos(phi), EZ);
  return E;
}

rad::ComplexVector3 rad::CircularCavity::GetModeEField(
    TVector3 pos, Mode_t modeType, double A, unsigned int n, unsigned int m,
    unsigned int l, bool state, double t) {
  double rho{pos.Perp()};
  double phi{atan2(pos.Y(), pos.X())};
  double z{pos.Z()};
  return GetModeEField(rho, phi, z, modeType, A, n, m, l, state, t);
}

rad::ComplexVector3 rad::CircularCavity::GetModalEField(
    TVector3 pos, Mode_t modeType, double A, unsigned int m, unsigned int n,
    unsigned int p, bool state) {
  // Get variables for calculation in terms of the formular
  const double z{pos.Z() + d / 2};
  const double phi{atan2(pos.Y(), pos.X())};
  double rho{pos.Perp()};

  // If we're outside the cavity then the field should be 0
  if (rho > a || z > d || z < 0) return ComplexVector3(0, 0, 0);

  std::complex<double> EPhi{0};
  std::complex<double> ERho{0};
  std::complex<double> EZ{0};

  const double beta{double(p) * TMath::Pi() / d};
  if (modeType == CircularCavity::kTE) {
    const double XmnPrime{GetBesselPrimeZero(m, n)};
    ERho = (-A * double(m) / rho) * sin(beta * z) *
           boost::math::cyl_bessel_j(m, XmnPrime * rho / a);
    EPhi = (A * XmnPrime / a) * sin(beta * z) *
           boost::math::cyl_bessel_j_prime(m, XmnPrime * rho / a);
    if (state) {
      ERho *= -sin(double(m) * phi);
      EPhi *= cos(double(m) * phi);
    } else {
      ERho *= cos(double(m) * phi);
      EPhi *= sin(double(m) * phi);
    }
  } else if (modeType == CircularCavity::kTM) {
    const double Xmn{boost::math::cyl_bessel_j_zero(double(m), n)};
    const double kc{Xmn / a};
    EZ = A * cos(beta * z) * boost::math::cyl_bessel_j(m, Xmn * rho / a);
    ERho = (-beta * Xmn / (a * kc * kc)) * sin(beta * z) *
           boost::math::cyl_bessel_j_prime(m, Xmn * rho / a);
    EPhi = (-A * beta * double(m) / (rho * kc * kc)) * sin(beta * z) *
           boost::math::cyl_bessel_j(m, Xmn * rho / a);
    if (state) {
      EZ *= cos(double(m) * phi);
      ERho *= cos(double(m) * phi);
      EPhi *= -sin(double(m) * phi);
    } else {
      EZ *= sin(double(m) * phi);
      ERho *= sin(double(m) * phi);
      EPhi *= cos(double(m) * phi);
    }
  } else {
    std::cout << "TEM modes not supported for circular cavities!\n";
  }

  // Generate final field vector
  ComplexVector3 E(ERho * cos(phi) - EPhi * sin(phi),
                   ERho * sin(phi) + EPhi * cos(phi), EZ);
  return E;
}

std::complex<double> rad::CircularCavity::GetModeNormalisation(Mode_t modeType,
                                                               unsigned int m,
                                                               unsigned int n,
                                                               unsigned int p,
                                                               bool state) {
  const unsigned int nIntPnts{50};
  const double dRho{a / nIntPnts};
  const double dPhi{TMath::TwoPi() / nIntPnts};
  const double dZ{d / nIntPnts};
  std::complex<double> volumeIntegral{0};
  // Integrate by looping over rho, phi and z
  for (unsigned int iRho{0}; iRho < nIntPnts; iRho++) {
    const double rho{dRho / 2 + double(iRho) * dRho};
    const double dV{rho * dRho * dPhi * dZ};
    for (unsigned int iPhi{0}; iPhi < nIntPnts; iPhi) {
      const double phi{dPhi / 2 + double(iPhi) * dPhi};
      for (unsigned int iZ{0}; iZ < nIntPnts; iZ++) {
        const double z{dZ / 2 + double(iZ) * dZ};
        TVector3 pos(rho * cos(phi), rho * sin(phi), z);

        // Now calculate the electric field at this point
        ComplexVector3 eField{GetModalEField(pos, modeType, 1, m, n, p, state)};
        volumeIntegral += eField.Dot(eField) * dV;
      }
    }
  }
  return 1.0 / sqrt(volumeIntegral);
}

rad::ComplexVector3 rad::CircularCavity::GetMaxEField(
    TVector3 pos, Mode_t modeType, double A, unsigned int m, unsigned int n,
    unsigned int p, bool state) {
  if (modeType == kTE) {
    // Calculate the resonant frequency
    const double fRes{GetResonantModeF(CircularCavity::kTE, m, n, p)};
    // Now calculate the time the field is at a maximum
    const double tMax{1 / (4 * fRes)};
    return GetModeEField(pos, kTE, A, m, n, p, state, tMax);
  } else if (modeType == kTM) {
    return GetModeEField(pos, kTM, A, m, n, p, state, 0);
  } else {
    std::cout << "TEM modes not supported for circular cavities! Returning "
                 "null vector\n";
    return ComplexVector3(0, 0, 0);
  }
}

rad::ComplexVector3 rad::CircularCavity::GetModeHField(
    double rho, double phi, double z, Mode_t modeType, double A, unsigned int m,
    unsigned int n, unsigned int p, bool state, double t) {
  const double f{GetResonantModeF(modeType, m, n, p)};
  const double omega{f * TMath::TwoPi()};
  const double z0{d / 2};

  // Check we're inside the cavity, return 0 if not
  if (rho > a || abs(z) > z0) return ComplexVector3(0, 0, 0);

  std::complex<double> HPhi{0};
  std::complex<double> HRho{0};
  std::complex<double> HZ{0};

  // Check we are inside the cavity
  if (modeType == CircularCavity::kTE) {
    const double prefactor{A / (TMath::C() * MU0)};
    const double XmnPrime{GetBesselPrimeZero(m, n)};
    const double factor1{cos(double(p) * TMath::PiOver2() * (z / z0 + 1)) *
                         pow(a / XmnPrime, 2) * double(p) * TMath::PiOver2() /
                         z0};
    if (state) {
      HZ = boost::math::cyl_bessel_j(m, XmnPrime * rho / a) *
           sin(double(p) * TMath::PiOver2() * (z / z0 + 1)) *
           cos(omega * t - double(m) * phi);
      HZ *= prefactor;
      HRho = (XmnPrime / a) *
             boost::math::cyl_bessel_j_prime(m, XmnPrime * rho / a) *
             cos(omega * t - double(m) * phi);
      HRho *= prefactor * factor1;
      HPhi = (double(m) / rho) *
             boost::math::cyl_bessel_j(m, XmnPrime * rho / a) *
             sin(omega * t - double(m) * phi);
      HPhi *= prefactor * factor1;
    } else {
      HZ = boost::math::cyl_bessel_j(m, XmnPrime * rho / a) *
           sin(double(p) * TMath::PiOver2() * (z / z0 + 1)) *
           cos(omega * t + double(m) * phi);
      HZ *= prefactor;
      HRho = (XmnPrime / a) *
             boost::math::cyl_bessel_j_prime(m, XmnPrime * rho / a) *
             cos(omega * t + double(m) * phi);
      HRho *= cos(double(p) * TMath::PiOver2() * (z / z0 + 1)) *
              pow(a / XmnPrime, 2) * double(p) * TMath::PiOver2() / z0;
      HRho *= prefactor;
      HPhi = (-double(m) / rho) *
             boost::math::cyl_bessel_j(m, XmnPrime * rho / a) *
             sin(omega * t + double(m) * phi);
      HPhi *= prefactor * factor1;
    }
  } else if (modeType == CircularCavity::kTM) {
    const double Xmn{boost::math::cyl_bessel_j_zero(double(m), n)};
    const double prefactor{(A * omega / (pow(TMath::C(), 2) * MU0)) *
                           pow(a / Xmn, 2) *
                           cos(double(p) * TMath::PiOver2() * (z / z0 + 1))};
    // Choose polarisation state
    if (state) {
      HRho = (double(m) / rho) * boost::math::cyl_bessel_j(m, Xmn * rho / a) *
             cos(omega * t - double(m) * phi);
      HRho *= prefactor;
      HPhi = (Xmn / a) * boost::math::cyl_bessel_j_prime(m, Xmn * rho / a) *
             sin(omega * t - double(m) * phi);
      HPhi *= prefactor;
    } else {
      HRho = (-double(m) / rho) * boost::math::cyl_bessel_j(m, Xmn * rho / a) *
             cos(omega * t + double(m) * phi);
      HRho *= prefactor;
      HPhi = (Xmn / a) * boost::math::cyl_bessel_j_prime(m, Xmn * rho / a) *
             sin(omega * t + double(m) * phi);
      HPhi *= prefactor;
    }
  }

  // Generate final field vector
  ComplexVector3 H(HRho * cos(phi) - HPhi * sin(phi),
                   HRho * sin(phi) + HPhi * cos(phi), HZ);
  return H;
}

rad::ComplexVector3 rad::CircularCavity::GetModeHField(
    TVector3 pos, Mode_t modeType, double A, unsigned int m, unsigned int n,
    unsigned int p, bool state, double t) {
  double rho{pos.Perp()};
  double phi{atan2(pos.Y(), pos.X())};
  double z{pos.Z()};
  return GetModeHField(rho, phi, z, modeType, A, m, n, p, state, t);
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