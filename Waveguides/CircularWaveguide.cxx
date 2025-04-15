/// CircularWaveguide.cxx

#include "Waveguides/CircularWaveguide.h"

#include <boost/math/special_functions/bessel.hpp>
#include <boost/math/special_functions/bessel_prime.hpp>
#include <cmath>
#include <complex>
#include <iostream>

#include "BasicFunctions/BasicFunctions.h"
#include "TMath.h"
#include "TVector3.h"

rad::CircularWaveguide::CircularWaveguide(double radius, double length)
    : a(radius), d(radius) {}

TVector3 rad::CircularWaveguide::GetModeEField(TVector3 pos, WaveguideMode mode,
                                               double A, double omega,
                                               bool state) {
  if (omega < GetCutoffFrequency(mode)) {
    std::cout
        << "Mode is below cutoff frequency and will not propagate. Returning 0."
        << std::endl;
    return TVector3(0, 0, 0);
  }
  unsigned int n{mode.GetModeIndex1()};
  unsigned int m{mode.GetModeIndex2()};

  double rho{pos.Perp()};
  if (rho > a) return TVector3(0, 0, 0);

  double phi{pos.Phi()};
  double ERho{0};
  double EPhi{0};
  double EZ{0};
  if (mode.GetModeType() == ModeType::kTE) {
    // We have a TE mode
    double XnmPrime{GetBesselPrimeZero(n, m)};
    ERho = (-A * double(n) / rho) *
           boost::math::cyl_bessel_j(n, XnmPrime * rho / a);
    EPhi = (A * XnmPrime / a) *
           boost::math::cyl_bessel_j_prime(n, XnmPrime * rho / a);
    if (state) {
      ERho *= cos(double(n) * phi);
      EPhi *= sin(double(n) * phi);
    } else {
      ERho *= -sin(double(n) * phi);
      EPhi *= cos(double(n) * phi);
    }
  } else if (mode.GetModeType() == ModeType::kTM) {
    // We have a TM mode
    double Xnm{boost::math::cyl_bessel_j_zero(double(n), m)};
    double k_c{Xnm / a};
    double betaSq{pow(omega / TMath::C(), 2) - k_c * k_c};
    double beta{sqrt(betaSq)};
    EZ = A * boost::math::cyl_bessel_j(n, k_c * rho);
    ERho = (-A * beta * Xnm / (a * k_c * k_c)) *
           boost::math::cyl_bessel_j_prime(n, k_c * rho);
    EPhi = (-A * beta * double(n) / (k_c * k_c * rho)) *
           boost::math::cyl_bessel_j(n, k_c * rho);
    if (state) {
      EZ *= sin(double(n) * phi);
      ERho *= sin(double(n) * phi);
      EPhi *= cos(double(n) * phi);
    } else {
      EZ *= cos(double(n) * phi);
      ERho *= cos(double(n) * phi);
      EPhi *= -sin(double(n) * phi);
    }
  } else {
    // This is a TEM mode that which is not supported by this waveguide
    std::cout << "TEM modes are not supported by circular waveguides."
              << std::endl;
  }
  TVector3 eField{ERho * cos(phi) - EPhi * sin(phi),
                  ERho * sin(phi) + EPhi * cos(phi), EZ};
  return eField;
}

TVector3 rad::CircularWaveguide::GetModeHField(TVector3 pos, WaveguideMode mode,
                                               double A, double omega,
                                               bool state) {
  if (omega < GetCutoffFrequency(mode)) {
    std::cout
        << "Mode is below cutoff frequency and will not propagate. Returning 0."
        << std::endl;
    return TVector3(0, 0, 0);
  }
  if (pos.Perp() > a) return TVector3(0, 0, 0);
  unsigned int n{mode.GetModeIndex1()};
  unsigned int m{mode.GetModeIndex2()};
  double rho{pos.Perp()};
  double phi{pos.Phi()};
  double HRho{0};
  double HPhi{0};
  double HZ{0};
  if (mode.GetModeType() == ModeType::kTE) {
    // We have a TE mode
    double XnmPrime{GetBesselPrimeZero(n, m)};
    double k_c{XnmPrime / a};
    double betaSq{pow(omega / TMath::C(), 2) - k_c * k_c};
    double beta{sqrt(betaSq)};
    HZ = -A * k_c * k_c * boost::math::cyl_bessel_j(n, k_c * rho) /
         (omega * MU0);
    HRho = -A * beta * k_c * boost::math::cyl_bessel_j_prime(n, k_c * rho) /
           (omega * MU0);
    HPhi = -A * beta * double(n) * boost::math::cyl_bessel_j(n, k_c * rho) /
           (rho * omega * MU0);
    if (state) {
      HZ *= sin(double(n) * phi);
      HRho *= sin(double(n) * phi);
      HPhi *= cos(double(n) * phi);
    } else {
      HZ *= cos(double(n) * phi);
      HRho *= cos(double(n) * phi);
      HPhi *= -sin(double(n) * phi);
    }
  } else if (mode.GetModeType() == ModeType::kTM) {
    // We have a TM mode
    double Xnm{boost::math::cyl_bessel_j_zero(double(n), m)};
    double k_c{Xnm / a};
    double betaSq{pow(omega / TMath::C(), 2) - k_c * k_c};
    double beta{sqrt(betaSq)};
    HRho = A * EPSILON0 * double(n) * boost::math::cyl_bessel_j(n, k_c * rho) /
           (rho * MU0);
    HPhi = -A * EPSILON0 * k_c * boost::math::cyl_bessel_j_prime(n, k_c * rho) /
           MU0;
    if (state) {
      HRho *= cos(double(n) * phi);
      HPhi *= sin(double(n) * phi);
    } else {
      HRho *= -sin(double(n) * phi);
      HPhi *= cos(double(n) * phi);
    }
  } else {
    // This is a TEM mode that which is not supported by this waveguide
    std::cout << "TEM modes are not supported by circular waveguides."
              << std::endl;
  }

  return TVector3(HRho * cos(phi) - HPhi * sin(phi),
                  HRho * sin(phi) + HPhi * cos(phi), HZ);
}

double rad::CircularWaveguide::GetCutoffFrequency(WaveguideMode mode) {
  ModeType mt{mode.GetModeType()};
  unsigned int n{mode.GetModeIndex1()};
  unsigned int m{mode.GetModeIndex2()};
  // Invalid mode number
  if (m == 0) {
    return -1;
  } else {
    if (mt == ModeType::kTE) {
      double pnmPrime{GetBesselPrimeZero(n, m)};
      double k_c{pnmPrime / a};
      double f_c{k_c * TMath::C() / (2 * TMath::Pi())};
      return f_c;
    } else if (mt == ModeType::kTM) {
      double pnm{boost::math::cyl_bessel_j_zero(double(n), m)};
      double k_c{pnm / a};
      double f_c{k_c * TMath::C() / (2 * TMath::Pi())};
      return f_c;
    } else {
      return -1;
    }
  }
}

double rad::CircularWaveguide::GetCutoffWavenumber(WaveguideMode mode) {
  ModeType mt{mode.GetModeType()};
  unsigned int n{mode.GetModeIndex1()};
  unsigned int m{mode.GetModeIndex2()};
  // Invalid mode number
  if (m == 0) {
    return -1;
  } else {
    if (mt == ModeType::kTE) {
      double pnmPrime{GetBesselPrimeZero(n, m)};
      double k_c{pnmPrime / a};
      return k_c;
    } else if (mt == ModeType::kTM) {
      double pnm{boost::math::cyl_bessel_j_zero(double(n), m)};
      double k_c{pnm / a};
      return k_c;
    } else {
      // This is a TEM mode that which is not supported by this waveguide
      std::cout << "TEM modes are not supported by circular waveguides."
                << std::endl;
      return -1;
    }
  }
}

double rad::CircularWaveguide::GetEFieldIntegral(WaveguideMode mode,
                                                 double omega, double A,
                                                 int nSurfPnts, bool state) {
  double integral{0};
  const double dRho{a / double(nSurfPnts)};
  const double dPhi{TMath::TwoPi() / double(nSurfPnts)};

  for (unsigned int iRho{0}; iRho < nSurfPnts; iRho++) {
    double thisRho{a / (2.0 * double(nSurfPnts)) +
                   a * double(iRho) / double(nSurfPnts)};
    const double area{thisRho * dRho * dPhi};
    for (unsigned int iPhi{0}; iPhi < nSurfPnts; iPhi++) {
      double thisPhi{TMath::TwoPi() / (2.0 * double(nSurfPnts)) +
                     TMath::TwoPi() * double(iPhi) / double(nSurfPnts)};

      TVector3 surfacePos{thisRho * cos(thisPhi), thisRho * sin(thisPhi), 0.0};
      TVector3 eTrans{GetModeEField(surfacePos, mode, A, omega, state)};
      eTrans.SetZ(0);
      integral += eTrans.Dot(eTrans) * area;
    }
  }

  if (integral == 0) {
    return 0;
  } else {
    return integral;
  }
}

double rad::CircularWaveguide::GetHFieldIntegral(WaveguideMode mode,
                                                 double omega, double A,
                                                 double B, int nSurfPnts) {
  const double dRho{a / double(nSurfPnts)};
  const double dPhi{TMath::TwoPi() / double(nSurfPnts)};
  double integral{0};
  for (unsigned int iRho{0}; iRho < nSurfPnts; iRho++) {
    double thisRho{a / (2.0 * double(nSurfPnts)) +
                   a * double(iRho) / double(nSurfPnts)};
    const double area{thisRho * dRho * dPhi};
    for (unsigned int iPhi{0}; iPhi < nSurfPnts; iPhi++) {
      double thisPhi{TMath::TwoPi() / (2.0 * double(nSurfPnts)) +
                     TMath::TwoPi() * double(iPhi) / double(nSurfPnts)};

      TVector3 surfacePos{thisRho * cos(thisPhi), thisRho * sin(thisPhi), 0.0};
      ComplexVector3 hTrans{GetModeHField(surfacePos, mode, A, omega, true)};
      hTrans.SetZ(std::complex<double>{0, 0});
      integral += (hTrans.Dot(hTrans)).real() * area;
    }
  }

  return integral;
}

double rad::CircularWaveguide::GetFieldAmp(WaveguideMode mode, double omega,
                                           TVector3 ePos, TVector3 eVel,
                                           double normA, bool state,
                                           bool isPositive) {
  TVector3 j{-TMath::Qe() * eVel};
  TVector3 eTrans{GetModeEField(ePos, mode, normA, omega, state)};
  eTrans.SetZ(0);
  TVector3 eAxial{GetModeEField(ePos, mode, normA, omega, state)};
  eAxial.SetX(0);
  eAxial.SetY(0);
  TVector3 eField(0, 0, 0);
  if (isPositive) {
    eField = eTrans - eAxial;
  } else {
    eField = eTrans + eAxial;
  }

  double A{eField.Dot(j)};
  return A;
}

void rad::CircularWaveguide::CalculatePn(WaveguideMode mode, double omega,
                                         unsigned int nSurfPnts) {
  double sum{0};
  // Calculate element of area
  for (unsigned int iRho{0}; iRho < nSurfPnts; iRho++) {
    double thisRho{GetInnerRadius() / (2 * double(nSurfPnts)) +
                   GetInnerRadius() * double(iRho) / double(nSurfPnts)};
    const double elArea{thisRho * (GetInnerRadius() / double(nSurfPnts)) *
                        (TMath::TwoPi() / double(nSurfPnts))};
    for (unsigned int iPhi{0}; iPhi < nSurfPnts; iPhi++) {
      double thisPhi{TMath::TwoPi() / (2 * double(nSurfPnts)) +
                     TMath::TwoPi() * double(iPhi) / double(nSurfPnts)};
      TVector3 surfacePos(thisRho * cos(thisPhi), thisRho * sin(thisPhi),
                          -GetLength() / 2);

      // Get transverse E and H components
      TVector3 eTrans{GetModeEField(surfacePos, mode, 1, omega, true)};
      eTrans.SetZ(0);
      TVector3 hTrans{GetModeHField(surfacePos, mode, 1, omega, true)};
      hTrans.SetZ(0);
      sum += (eTrans.Cross(hTrans)).Dot(TVector3(0, 0, 1)) * elArea;
    }
  }
  SetPn(sum * 2);
}