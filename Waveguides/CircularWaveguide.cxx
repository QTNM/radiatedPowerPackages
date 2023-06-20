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

rad::ComplexVector3 rad::CircularWaveguide::GetModeEFieldComplex(
    Mode_t modeType, int n, int m, TVector3 pos, double omega, double A,
    double B) {
  double rho{pos.Perp()};
  double phi{pos.Phi()};
  double z{pos.Z() + d / 2.0};
  std::complex<double> i{0.0, 1.0};

  if (modeType == kTE) {
    // We have a TE mode
    double pnmPrime{GetBesselPrimeZero(n, m)};
    double k_c{pnmPrime / a};
    std::complex<double> betaSq{pow(omega / TMath::C(), 2) - k_c * k_c};
    std::complex<double> beta{sqrt(betaSq)};

    std::complex<double> Ez{0.0, 0.0};
    std::complex<double> Erho{
        (-1.0 * omega * MU0 * double(n) / (k_c * k_c * rho)) *
        (A * cos(double(n) * phi) - B * sin(double(n) * phi)) *
        boost::math::cyl_bessel_j(n, k_c * rho)};
    Erho *= i * exp(-1.0 * i * beta * z);
    std::complex<double> Ephi{
        (omega * MU0 / k_c) *
        (A * sin(double(n) * phi) + B * cos(double(n) * phi)) *
        boost::math::cyl_bessel_j_prime(n, k_c * rho)};
    Ephi *= i * exp(-1.0 * i * beta * z);

    ComplexVector3 eField{Erho * cos(phi) - Ephi * sin(phi),
                          Erho * sin(phi) + Ephi * cos(phi), Ez};
    return eField;
  } else if (modeType == kTM) {
    // We have a TM mode
    double pnm{boost::math::cyl_bessel_j_zero(double(n), m)};
    double k_c{pnm / a};
    std::complex<double> betaSq{pow(omega / TMath::C(), 2) - k_c * k_c};
    std::complex<double> beta{sqrt(betaSq)};

    std::complex<double> Ez{
        (A * sin(double(n) * phi) + B * cos(double(n) * phi)) *
            boost::math::cyl_bessel_j(n, k_c * rho),
        0.0};
    Ez *= exp(-1.0 * i * beta * z);
    std::complex<double> Erho{
        (-1.0 * beta / k_c) *
        (A * sin(double(n) * phi) + B * cos(double(n) * phi)) *
        boost::math::cyl_bessel_j_prime(n, k_c * rho)};
    Erho *= i * exp(-1.0 * i * beta * z);
    std::complex<double> Ephi{
        (-1.0 * beta * double(n) / (k_c * k_c * rho)) *
        (A * cos(double(n) * phi) - B * sin(double(n) * phi)) *
        boost::math::cyl_bessel_j(n, k_c * rho)};
    Ephi *= i * exp(-1.0 * i * beta * z);

    ComplexVector3 eField{Erho * cos(phi) - Ephi * sin(phi),
                          Erho * sin(phi) + Ephi * cos(phi), Ez};
    return eField;
  } else {
    // This is a TEM mode that which is not supported by this waveguide
    std::cout << "TEM modes are not supported by circular waveguides."
              << std::endl;
    ComplexVector3 eField{std::complex<double>{0}, std::complex<double>{0},
                          std::complex<double>{0}};
    return eField;
  }
}

TVector3 rad::CircularWaveguide::GetModeEField(Mode_t modeType, int n, int m,
                                               TVector3 pos, double omega,
                                               double A, double B) {
  ComplexVector3 field{GetModeEFieldComplex(modeType, n, m, pos, omega, A, B)};
  return field.Real();
}

rad::ComplexVector3 rad::CircularWaveguide::GetModalEField(Mode_t modeType,
                                                           int n, int m,
                                                           TVector3 pos,
                                                           double omega,
                                                           double A, double B) {
  double rho{pos.Perp()};
  double phi{pos.Phi()};
  double z{pos.Z() + d / 2.0};
  std::complex<double> i{0.0, 1.0};

  if (modeType == kTE) {
    // We have a TE mode
    double pnmPrime{GetBesselPrimeZero(n, m)};
    double k_c{pnmPrime / a};
    std::complex<double> betaSq{pow(omega / TMath::C(), 2) - k_c * k_c};
    std::complex<double> beta{sqrt(betaSq)};
    std::complex<double> Ez{0.0, 0.0};
    std::complex<double> Erho{
        (omega * MU0 * double(n) / (k_c * k_c * rho)) *
        (A * cos(double(n) * phi) - B * sin(double(n) * phi)) *
        boost::math::cyl_bessel_j(n, k_c * rho)};
    std::complex<double> Ephi{
        (-1.0 * omega * MU0 / k_c) *
        (A * sin(double(n) * phi) + B * cos(double(n) * phi)) *
        boost::math::cyl_bessel_j_prime(n, k_c * rho)};

    ComplexVector3 eField{Erho * cos(phi) - Ephi * sin(phi),
                          Erho * sin(phi) + Ephi * cos(phi), Ez};
    return eField;
  } else if (modeType == kTM) {
    // We have a TM mode
    double pnm{boost::math::cyl_bessel_j_zero(double(n), m)};
    double k_c{pnm / a};
    std::complex<double> betaSq{pow(omega / TMath::C(), 2) - k_c * k_c};
    std::complex<double> beta{sqrt(betaSq)};
    std::complex<double> Ez{
        i * (A * sin(double(n) * phi) + B * cos(double(n) * phi)) *
        boost::math::cyl_bessel_j(n, k_c * rho)};
    std::complex<double> Erho{
        (beta / k_c) * (A * sin(double(n) * phi) + B * cos(double(n) * phi)) *
        boost::math::cyl_bessel_j_prime(n, k_c * rho)};
    std::complex<double> Ephi{
        (beta * double(n) / (k_c * k_c * rho)) *
        (A * cos(double(n) * phi) - B * sin(double(n) * phi)) *
        boost::math::cyl_bessel_j(n, k_c * rho)};

    ComplexVector3 eField{Erho * cos(phi) - Ephi * sin(phi),
                          Erho * sin(phi) + Ephi * cos(phi), Ez};
    return eField;
  } else {
    // This is a TEM mode that which is not supported by this waveguide
    std::cout << "TEM modes are not supported by circular waveguides."
              << std::endl;
    ComplexVector3 eField{std::complex<double>{0, 0},
                          std::complex<double>{0, 0},
                          std::complex<double>{0, 0}};
    return eField;
  }
}

rad::ComplexVector3 rad::CircularWaveguide::GetModeHFieldComplex(
    Mode_t modeType, int n, int m, TVector3 pos, double omega, double A,
    double B) {
  double rho{pos.Perp()};
  double phi{pos.Phi()};
  double z{pos.Z() + d / 2.0};
  std::complex<double> i{0.0, 1.0};

  if (modeType == kTE) {
    // We have a TE mode
    double pnmPrime{GetBesselPrimeZero(n, m)};
    double k_c{pnmPrime / a};
    std::complex<double> betaSq{pow(omega / TMath::C(), 2) - k_c * k_c};
    std::complex<double> beta{sqrt(betaSq)};

    std::complex<double> Hz{
        (A * sin(double(n) * phi) + B * cos(double(n) * phi)) *
            boost::math::cyl_bessel_j(n, k_c * rho),
        0.0};
    Hz *= exp(-1.0 * i * beta * z);
    std::complex<double> Hrho{
        (-1.0 * beta / k_c) *
        (A * sin(double(n) * phi) + B * cos(double(n) * phi)) *
        boost::math::cyl_bessel_j_prime(n, k_c * rho)};
    Hrho *= i * exp(-1.0 * i * beta * z);
    std::complex<double> Hphi{
        (-1.0 * beta * double(n) / (k_c * k_c * rho)) *
        (A * cos(double(n) * phi) - B * sin(double(n) * phi)) *
        boost::math::cyl_bessel_j(n, k_c * rho)};
    Hphi *= i * exp(-1.0 * i * beta * z);

    ComplexVector3 hField{Hrho * cos(phi) - Hphi * sin(phi),
                          Hrho * sin(phi) + Hphi * cos(phi), Hz};
    return hField;
  } else if (modeType == kTM) {
    // We have a TM mode
    double pnm{boost::math::cyl_bessel_j_zero(double(n), m)};
    double k_c{pnm / a};
    std::complex<double> betaSq{pow(omega / TMath::C(), 2) - k_c * k_c};
    std::complex<double> beta{sqrt(betaSq)};

    std::complex<double> Hz{0.0, 0.0};
    std::complex<double> Hrho{
        (omega * EPSILON0 * double(n) / (k_c * k_c * rho)) *
        (A * cos(double(n) * phi) - B * sin(double(n) * phi)) *
        boost::math::cyl_bessel_j(n, k_c * rho)};
    Hrho *= i * exp(-1.0 * i * beta * z);
    std::complex<double> Hphi{
        (-1.0 * omega * EPSILON0 / k_c) *
        (A * sin(double(n) * phi) + B * cos(double(n) * phi)) *
        boost::math::cyl_bessel_j_prime(n, k_c * rho)};
    Hphi *= i * exp(-1.0 * i * beta * z);

    ComplexVector3 hField{Hrho * cos(phi) - Hphi * sin(phi),
                          Hrho * sin(phi) + Hphi * cos(phi), Hz};
    return hField;
  } else {
    // This is a TEM mode that which is not supported by this waveguide
    std::cout << "TEM modes are not supported by circular waveguides."
              << std::endl;
    ComplexVector3 eField{std::complex<double>{0}, std::complex<double>{0},
                          std::complex<double>{0}};
    return eField;
  }
}

rad::ComplexVector3 rad::CircularWaveguide::GetModalHField(Mode_t modeType,
                                                           int n, int m,
                                                           TVector3 pos,
                                                           double omega,
                                                           double A, double B) {
  double rho{pos.Perp()};
  double phi{pos.Phi()};
  double z{pos.Z() + d / 2.0};
  std::complex<double> i{0.0, 1.0};

  if (modeType == kTE) {
    // We have a TE mode
    double pnmPrime{GetBesselPrimeZero(n, m)};
    double k_c{pnmPrime / a};
    std::complex<double> betaSq{pow(omega / TMath::C(), 2) - k_c * k_c};
    std::complex<double> beta{sqrt(betaSq)};

    std::complex<double> Hrho{
        (beta / k_c) * (A * sin(double(n) * phi) + B * cos(double(n) * phi)) *
        boost::math::cyl_bessel_j_prime(n, k_c * rho)};
    std::complex<double> Hphi{
        (beta * double(n) / (k_c * k_c * rho)) *
        (A * cos(double(n) * phi) - B * sin(double(n) * phi)) *
        boost::math::cyl_bessel_j(n, k_c * rho)};
    std::complex<double> Hz{
        i * (A * sin(double(n) * phi) + B * cos(double(n) * phi)) *
        boost::math::cyl_bessel_j(n, k_c * rho)};

    ComplexVector3 hField{Hrho * cos(phi) - Hphi * sin(phi),
                          Hrho * sin(phi) + Hphi * cos(phi), Hz};
    return hField;
  } else if (modeType == kTM) {
    // We have a TM mode
    double pnm{boost::math::cyl_bessel_j_zero(double(n), m)};
    double k_c{pnm / a};
    std::complex<double> betaSq{pow(omega / TMath::C(), 2) - k_c * k_c};
    std::complex<double> beta{sqrt(betaSq)};

    std::complex<double> Hrho{
        (-1.0 * omega * EPSILON0 * double(n) / (k_c * k_c * rho)) *
        (A * cos(double(n) * phi) - B * sin(double(n) * phi)) *
        boost::math::cyl_bessel_j(n, k_c * rho)};
    std::complex<double> Hphi{
        (omega * EPSILON0 / k_c) *
        (A * sin(double(n) * phi) + B * cos(double(n) * phi)) *
        boost::math::cyl_bessel_j_prime(n, k_c * rho)};
    std::complex<double> Hz{0.0, 0.0};

    ComplexVector3 hField{Hrho * cos(phi) - Hphi * sin(phi),
                          Hrho * sin(phi) + Hphi * cos(phi), Hz};
    return hField;
  } else {
    // This is a TEM mode that which is not supported by this waveguide
    std::cout << "TEM modes are not supported by circular waveguides."
              << std::endl;
    ComplexVector3 eField{std::complex<double>{0}, std::complex<double>{0},
                          std::complex<double>{0}};
    return eField;
  }
}

TVector3 rad::CircularWaveguide::GetModeHField(Mode_t modeType, int n, int m,
                                               TVector3 pos, double omega,
                                               double A, double B) {
  ComplexVector3 field{GetModeHFieldComplex(modeType, n, m, pos, omega, A, B)};
  return field.Real();
}

double rad::CircularWaveguide::GetCutoffFrequency(Mode_t modeType, int n,
                                                  int m) {
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
    return -1;
  }
}

double rad::CircularWaveguide::GetCutoffWavenumber(Mode_t modeType,
                                                   unsigned int n,
                                                   unsigned int m) {
  if (modeType == kTE) {
    double pnmPrime{GetBesselPrimeZero(n, m)};
    double k_c{pnmPrime / a};
    return k_c;
  } else if (modeType == kTM) {
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

double rad::CircularWaveguide::GetResonantModeFrequency(Mode_t modeType, int n,
                                                        int m, int l) {
  double premult{TMath::C() / (2 * TMath::C())};
  if (modeType == kTE) {
    double pnmPrime{GetBesselPrimeZero(n, m)};
    double freq{premult * sqrt(pow(pnmPrime / a, 2) +
                               pow(double(l) * TMath::Pi() / d, 2))};
    return freq;
  } else if (modeType == kTM) {
    double pnm{boost::math::cyl_bessel_j_zero(double(n), m)};
    double freq{premult *
                sqrt(pow(pnm / a, 2) + pow(double(l) * TMath::Pi() / d, 2))};
    return freq;
  } else {
    return -1;
  }
}

double rad::CircularWaveguide::GetEFieldNormalisation(Mode_t modeType, int n,
                                                      int m, double omega,
                                                      double A, double B,
                                                      int nSurfPnts) {
  double area{pow(a / double(nSurfPnts), 2)};
  double integral{0};
  for (int ix = 0; ix < nSurfPnts; ix++) {
    double thisx{-a / 2.0 + a / (2.0 * double(nSurfPnts)) +
                 a * double(ix) / double(nSurfPnts)};
    for (int iy = 0; iy < nSurfPnts; iy++) {
      double thisy{-a / 2.0 + a / (2.0 * double(nSurfPnts)) +
                   a * double(iy) / double(nSurfPnts)};

      if (sqrt(thisx * thisx + thisy * thisy) > a) continue;

      TVector3 surfacePos{thisx, thisy, 0.0};
      ComplexVector3 eTrans{
          GetModalEField(modeType, n, m, surfacePos, omega, A, B)};
      eTrans.SetZ(std::complex<double>{0, 0});
      integral += (eTrans.Dot(eTrans)).real() * area;
    }
  }
  std::cout << "integral = " << integral << std::endl;

  if (integral == 0) {
    return 0;
  } else {
    return sqrt(1.0 / integral);
  }
}

double rad::CircularWaveguide::GetHFieldIntegral(Mode_t modeType, int n, int m,
                                                 double omega, double A,
                                                 double B, int nSurfPnts) {
  double area{pow(a / double(nSurfPnts), 2)};
  double integral{0};
  for (int ix = 0; ix < nSurfPnts; ix++) {
    double thisx{-a / 2.0 + a / (2.0 * double(nSurfPnts)) +
                 a * double(ix) / double(nSurfPnts)};
    for (int iy = 0; iy < nSurfPnts; iy++) {
      double thisy{-a / 2.0 + a / (2.0 * double(nSurfPnts)) +
                   a * double(iy) / double(nSurfPnts)};

      if (sqrt(thisx * thisx + thisy * thisy) > a) continue;

      TVector3 surfacePos{thisx, thisy, 0.0};
      ComplexVector3 hTrans{
          GetModalHField(modeType, n, m, surfacePos, omega, A, B)};
      hTrans.SetZ(std::complex<double>{0, 0});
      integral += (hTrans.Dot(hTrans)).real() * area;
    }
  }

  return integral;
}

std::complex<double> rad::CircularWaveguide::GetPositiveFieldAmp(
    Mode_t modeType, unsigned int n, unsigned int m, double omega,
    TVector3 ePos, TVector3 eVel, double normA, double normB) {
  double waveImp{GetModeImpedance(modeType, n, m, omega)};
  TVector3 j{-TMath::Qe() * eVel};
  ComplexVector3 jComplex{j};

  ComplexVector3 eTrans{
      GetModalEField(modeType, n, m, ePos, omega, normA, normB)};
  eTrans.SetZ(std::complex<double>{0.0, 0.0});
  ComplexVector3 eAxial{
      GetModalEField(modeType, n, m, ePos, omega, normA, normB)};
  eAxial.SetX(std::complex<double>{0.0, 0.0});
  eAxial.SetY(std::complex<double>{0.0, 0.0});
  ComplexVector3 subVec{eTrans - eAxial};

  std::complex<double> APlus = (subVec).Dot(jComplex) * (-waveImp / 2.0);
  return APlus;
}

void rad::CircularWaveguide::CalculatePn(Mode_t modeType, unsigned int n,
                                         unsigned int m, double omega,
                                         unsigned int nSurfPnts, double A,
                                         double B) {
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
      TVector3 sPos1(thisRho * cos(thisPhi), thisRho * sin(thisPhi),
                     -GetLength() / 2);
      TVector3 sPos2(thisRho * cos(thisPhi), thisRho * sin(thisPhi),
                     GetLength() / 2);

      // Get transverse E and H components
      ComplexVector3 eTrans1{
          GetModalEField(modeType, n, m, sPos1, omega, A, B)};
      eTrans1.SetZ(0);
      ComplexVector3 hTrans1{
          GetModalHField(modeType, n, m, sPos1, omega, A, B)};
      hTrans1.SetZ(0);
      ComplexVector3 eTrans2{
          GetModalEField(modeType, n, m, sPos2, omega, A, B)};
      eTrans2.SetZ(0);
      ComplexVector3 hTrans2{
          GetModalHField(modeType, n, m, sPos2, omega, A, B)};
      hTrans2.SetZ(0);
      sum += (eTrans1.Cross(hTrans1)).Z().real() * elArea;
      sum += (eTrans2.Cross(hTrans2)).Z().real() * elArea;
    }
  }
  SetPn(sum * 2);
}