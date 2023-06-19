/// RectangularWaveguide.cxx

#include "Waveguides/RectangularWaveguide.h"

#include <cmath>
#include <complex>
#include <iostream>

#include "BasicFunctions/Constants.h"
#include "TMath.h"
#include "TVector3.h"

rad::RectangularWaveguide::RectangularWaveguide(double longSide,
                                                double shortSide,
                                                double length) {
  d = length;
  if (longSide >= shortSide) {
    a = longSide;
    b = shortSide;
  } else {
    a = shortSide;
    b = longSide;
  }
}

double rad::RectangularWaveguide::GetCutoffWavenumber(Mode_t modeType,
                                                      unsigned int m,
                                                      unsigned int n) {
  double k_c{sqrt(pow(double(m) * TMath::Pi() / a, 2) +
                  pow(double(n) * TMath::Pi() / b, 2))};
  return k_c;
}

rad::ComplexVector3 rad::RectangularWaveguide::GetModeEFieldComplex(
    Mode_t modeType, int m, int n, TVector3 pos, double omega, double A,
    double B) {
  double x{pos.X() + a / 2.0};
  double y{pos.Y() + b / 2.0};
  double z{pos.Z() + d / 2.0};
  std::complex<double> i{0.0, 1.0};
  double k_c{GetCutoffWavenumber(modeType, m, n)};
  std::complex<double> betaSq{pow(omega / TMath::C(), 2) - k_c * k_c};
  std::complex<double> beta{sqrt(betaSq)};

  if (modeType == kTE) {
    // We have a TE mode
    std::complex<double> Ex{
        (i * omega * MU0 * double(n) * TMath::Pi() / (k_c * k_c * b)) * A *
        cos(double(m) * TMath::Pi() * x / a) *
        sin(double(n) * TMath::Pi() * y / b) * exp(-i * beta * z)};
    std::complex<double> Ey{
        (-1.0 * i * omega * MU0 * double(m) * TMath::Pi() / (k_c * k_c * a)) *
        A * sin(double(m) * TMath::Pi() * x / a) *
        cos(double(n) * TMath::Pi() * y / b) * exp(-i * beta * z)};
    std::complex<double> Ez{0.0, 0.0};
    ComplexVector3 eField{Ex, Ey, Ez};
    return eField;
  } else if (modeType == kTM) {
    // We have a TM mode
    std::complex<double> Ex{
        (-1.0 * i * beta * double(m) * TMath::Pi() / (k_c * k_c * a)) * A *
        cos(double(m) * TMath::Pi() * x / a) *
        sin(double(n) * TMath::Pi() * y / b) * exp(-i * beta * z)};
    std::complex<double> Ey{
        (-1.0 * i * beta * double(n) * TMath::Pi() / (k_c * k_c * b)) * A *
        sin(double(m) * TMath::Pi() * x / a) *
        cos(double(n) * TMath::Pi() * y / b) * exp(-i * beta * z)};
    std::complex<double> Ez{A * sin(double(m) * TMath::Pi() * x / a) *
                            sin(double(n) * TMath::Pi() * y / b) *
                            exp(-i * beta * z)};
    ComplexVector3 eField{Ex, Ey, Ez};
    return eField;
  } else {
    // We have a TEM mode
    std::cout << "TEM modes are not supported by circular waveguides."
              << std::endl;
    ComplexVector3 eField{std::complex<double>{0}, std::complex<double>{0},
                          std::complex<double>{0}};
    return eField;
  }
}

TVector3 rad::RectangularWaveguide::GetModeEField(Mode_t modeType, int m, int n,
                                                  TVector3 pos, double omega,
                                                  double A, double B) {
  ComplexVector3 field{GetModeEFieldComplex(modeType, m, n, pos, omega, A, B)};
  return field.Real();
}

TVector3 rad::RectangularWaveguide::GetModalEField(Mode_t modeType, int m,
                                                   int n, TVector3 pos,
                                                   double omega, double A) {
  double x{pos.X() + a / 2.0};
  double y{pos.Y() + b / 2.0};
  double z{pos.Z() + d / 2.0};
  double k_c{GetCutoffWavenumber(modeType, m, n)};
  double betaSq{pow(omega / TMath::C(), 2) - k_c * k_c};
  double beta{sqrt(betaSq)};

  if (modeType == kTE) {
    // We have a TE mode
    double Ex{(omega * MU0 * double(n) * TMath::Pi() / (k_c * k_c * b)) * A *
              cos(double(m) * TMath::Pi() * x / a) *
              sin(double(n) * TMath::Pi() * y / b)};
    double Ey{(-1.0 * omega * MU0 * double(m) * TMath::Pi() / (k_c * k_c * a)) *
              A * sin(double(m) * TMath::Pi() * x / a) *
              cos(double(n) * TMath::Pi() * y / b)};
    double Ez{0.0};
    TVector3 eField{Ex, Ey, Ez};
    return eField;
  } else if (modeType == kTM) {
    // We have a TM mode
    double Ex{(-1.0 * beta * double(m) * TMath::Pi() / (k_c * k_c * a)) * A *
              cos(double(m) * TMath::Pi() * x / a) *
              sin(double(n) * TMath::Pi() * y / b)};
    double Ey{(-1.0 * beta * double(n) * TMath::Pi() / (k_c * k_c * b)) * A *
              sin(double(m) * TMath::Pi() * x / a) *
              cos(double(n) * TMath::Pi() * y / b)};
    double Ez{A * sin(double(m) * TMath::Pi() * x / a) *
              sin(double(n) * TMath::Pi() * y / b)};
    TVector3 eField{Ex, Ey, Ez};
    return eField;
  } else {
    // We have a TEM mode
    std::cout << "TEM modes are not supported by rectangular waveguides."
              << std::endl;
    TVector3 eField{0, 0, 0};
    return eField;
  }
}

rad::ComplexVector3 rad::RectangularWaveguide::GetNormalisedEField(
    Mode_t modeType, int m, int n, TVector3 pos, double omega) {
  double x{pos.X() + a / 2.0};
  double y{pos.Y() + b / 2.0};
  double z{pos.Z() + d / 2.0};
  double k_c{GetCutoffWavenumber(modeType, m, n)};
  double betaSq{pow(omega / TMath::C(), 2) - k_c * k_c};
  double beta{sqrt(betaSq)};

  if (modeType == kTE) {
    std::complex<double> Ex{(-2.0 * TMath::Pi() * double(n)) /
                                (k_c * b * sqrt(a * b)) *
                                cos(double(m) * TMath::Pi() * x / a) *
                                sin(double(n) * TMath::Pi() * y / b),
                            0.0};
    std::complex<double> Ey{(-2.0 * TMath::Pi() * double(m)) /
                                (k_c * a * sqrt(a * b)) *
                                sin(double(m) * TMath::Pi() * x / a) *
                                cos(double(n) * TMath::Pi() * y / b),
                            0.0};
    std::complex<double> Ez{0, 0};
    ComplexVector3 eField{Ex, Ey, Ez};
    if (m == 0 || n == 0) eField *= (1.0 / sqrt(2.0));
    return eField;
  } else if (modeType == kTM) {
    std::complex<double> i{0, 1};
    std::complex<double> Ex{(2.0 * TMath::Pi() * double(m)) /
                                (k_c * a * sqrt(a * b)) *
                                cos(double(m) * TMath::Pi() * x / a) *
                                sin(double(n) * TMath::Pi() * y / b),
                            0};
    std::complex<double> Ey{(2.0 * TMath::Pi() * double(n)) /
                                (k_c * b * sqrt(a * b)) *
                                sin(double(m) * TMath::Pi() * x / a) *
                                cos(double(n) * TMath::Pi() * y / b),
                            0};
    std::complex<double> Ez{(-2.0 * i * k_c) / (beta * sqrt(a * b)) *
                            sin(double(m) * TMath::Pi() * x / a) *
                            sin(double(n) * TMath::Pi() * y / b)};
    ComplexVector3 eField{Ex, Ey, Ez};
    return eField;
  } else {
    // We have a TEM mode
    std::cout << "TEM modes are not supported by rectangular waveguides."
              << std::endl;
    ComplexVector3 eField{std::complex<double>{0, 0},
                          std::complex<double>{0, 0},
                          std::complex<double>{0, 0}};
    return eField;
  }
}

rad::ComplexVector3 rad::RectangularWaveguide::GetModeHFieldComplex(
    Mode_t modeType, int m, int n, TVector3 pos, double omega, double A,
    double B) {
  double x{pos.X() + a / 2.0};
  double y{pos.Y() + b / 2.0};
  double z{pos.Z() + d / 2.0};
  std::complex<double> i{0.0, 1.0};
  double k_c{GetCutoffWavenumber(modeType, m, n)};
  std::complex<double> betaSq{pow(omega / TMath::C(), 2) - k_c * k_c};
  std::complex<double> beta{sqrt(betaSq)};

  if (modeType == kTE) {
    // We have a TE mode
    std::complex<double> Hx{
        (i * beta * double(m) * TMath::Pi() / (k_c * k_c * a)) * A *
        sin(double(m) * TMath::Pi() * x / a) *
        cos(double(n) * TMath::Pi() * y / b) * exp(-i * beta * z)};
    std::complex<double> Hy{
        (i * beta * double(n) * TMath::Pi() / (k_c * k_c * b)) * A *
        cos(double(m) * TMath::Pi() * x / a) *
        sin(double(n) * TMath::Pi() * y / b) * exp(-i * beta * z)};
    std::complex<double> Hz{A * cos(double(m) * TMath::Pi() * x / a) *
                            cos(double(n) * TMath::Pi() * y / b) *
                            exp(-i * beta * z)};
    ComplexVector3 hField{Hx, Hy, Hz};
    return hField;
  } else if (modeType == kTM) {
    // We have a TM mode
    std::complex<double> Hx{
        (i * omega * EPSILON0 * double(n) * TMath::Pi() / (k_c * k_c * b)) * A *
        sin(double(m) * TMath::Pi() * x / a) *
        cos(double(n) * TMath::Pi() * y / b) * exp(-i * beta * z)};
    std::complex<double> Hy{(-1.0 * i * omega * EPSILON0 * double(m) *
                             TMath::Pi() / (k_c * k_c * a)) *
                            A * cos(double(m) * TMath::Pi() * x / a) *
                            sin(double(n) * TMath::Pi() * y / b) *
                            exp(-i * beta * z)};
    std::complex<double> Hz{0.0, 0.0};
    ComplexVector3 hField{Hx, Hy, Hz};
    return hField;
  } else {
    // We have a TEM mode
    std::cout << "TEM modes are not supported by circular waveguides."
              << std::endl;
    ComplexVector3 hField{std::complex<double>{0}, std::complex<double>{0},
                          std::complex<double>{0}};
    return hField;
  }
}

TVector3 rad::RectangularWaveguide::GetModeHField(Mode_t modeType, int m, int n,
                                                  TVector3 pos, double omega,
                                                  double A, double B) {
  ComplexVector3 field{GetModeHFieldComplex(modeType, m, n, pos, omega, A, B)};
  return field.Real();
}

TVector3 rad::RectangularWaveguide::GetModalHField(Mode_t modeType, int m,
                                                   int n, TVector3 pos,
                                                   double omega, double A) {
  double x{pos.X() + a / 2.0};
  double y{pos.Y() + b / 2.0};
  double z{pos.Z() + d / 2.0};
  double k_c{GetCutoffWavenumber(modeType, m, n)};
  double betaSq{pow(omega / TMath::C(), 2) - k_c * k_c};
  double beta{sqrt(betaSq)};

  if (modeType == kTE) {
    // We have a TE mode
    double Hx{(beta * double(m) * TMath::Pi() / (k_c * k_c * a)) * A *
              sin(double(m) * TMath::Pi() * x / a) *
              cos(double(n) * TMath::Pi() * y / b)};
    double Hy{(beta * double(n) * TMath::Pi() / (k_c * k_c * b)) * A *
              cos(double(m) * TMath::Pi() * x / a) *
              sin(double(n) * TMath::Pi() * y / b)};
    double Hz{A * cos(double(m) * TMath::Pi() * x / a) *
              cos(double(n) * TMath::Pi() * y / b)};
    TVector3 hField{Hx, Hy, Hz};
    return hField;
  } else if (modeType == kTM) {
    // We have a TM mode
    double Hx{(omega * EPSILON0 * double(n) * TMath::Pi() / (k_c * k_c * b)) *
              A * sin(double(m) * TMath::Pi() * x / a) *
              cos(double(n) * TMath::Pi() * y / b)};
    double Hy{
        (-1.0 * omega * EPSILON0 * double(m) * TMath::Pi() / (k_c * k_c * a)) *
        A * cos(double(m) * TMath::Pi() * x / a) *
        sin(double(n) * TMath::Pi() * y / b)};
    double Hz{0.0};
    TVector3 hField{Hx, Hy, Hz};
    return hField;
  } else {
    // We have a TEM mode
    std::cout << "TEM modes are not supported by rectangular waveguides."
              << std::endl;
    TVector3 hField{0, 0, 0};
    return hField;
  }
}

rad::ComplexVector3 rad::RectangularWaveguide::GetNormalisedHField(
    Mode_t modeType, int m, int n, TVector3 pos, double omega) {
  double x{pos.X() + a / 2.0};
  double y{pos.Y() + b / 2.0};
  double z{pos.Z() + d / 2.0};
  double k_c{GetCutoffWavenumber(modeType, m, n)};
  double betaSq{pow(omega / TMath::C(), 2) - k_c * k_c};
  double beta{sqrt(betaSq)};
  std::complex<double> i{0, 1};
  double waveImp{GetModeImpedance(modeType, m, n, omega)};
  TVector3 unitZ{0, 0, 1};
  ComplexVector3 unitZComplex{unitZ};

  if (modeType == kTE) {
    ComplexVector3 eFieldTrans{GetNormalisedEField(modeType, m, n, pos, omega)};
    eFieldTrans.SetZ(std::complex<double>{0, 0});
    ComplexVector3 hFieldTrans{unitZComplex.Cross(eFieldTrans) *
                               (1.0 / waveImp)};
    std::complex<double> Hz{(-2.0 * i * k_c) / (beta * waveImp * sqrt(a * b)) *
                            cos(double(m) * TMath::Pi() * x / a) *
                            cos(double(n) * TMath::Pi() * y / b)};
    ComplexVector3 hField{hFieldTrans.X(), hFieldTrans.Y(), Hz};
    return hField;
  } else if (modeType == kTM) {
    ComplexVector3 eFieldTrans{GetNormalisedEField(modeType, m, n, pos, omega)};
    eFieldTrans.SetZ(std::complex<double>{0, 0});
    ComplexVector3 hFieldTrans{unitZComplex.Cross(eFieldTrans) *
                               (1.0 / waveImp)};
    std::complex<double> Hz{0, 0};
    ComplexVector3 hField{hFieldTrans.X(), hFieldTrans.Y(), Hz};
    return hField;
  } else {
    // We have a TEM mode
    std::cout << "TEM modes are not supported by rectangular waveguides."
              << std::endl;
    ComplexVector3 hField{std::complex<double>{0, 0},
                          std::complex<double>{0, 0},
                          std::complex<double>{0, 0}};
    return hField;
  }
}

double rad::RectangularWaveguide::GetCutoffFrequency(Mode_t modeType, int m,
                                                     int n) {
  double k_c{GetCutoffWavenumber(modeType, m, n)};
  return k_c * TMath::C() / (2 * TMath::Pi());
}

double rad::RectangularWaveguide::GetResonantModeFrequency(Mode_t modeType,
                                                           int m, int n,
                                                           int l) {
  double freq{(TMath::C() / (2 * TMath::Pi())) *
              sqrt(pow(double(m) * TMath::Pi() / a, 2) +
                   pow(double(n) * TMath::Pi() / b, 2) +
                   pow(double(l) * TMath::Pi() / d, 2))};
  return freq;
}

std::complex<double> rad::RectangularWaveguide::GetPositiveFieldAmp(
    Mode_t modeType, unsigned int m, unsigned int n, double omega,
    TVector3 ePos, TVector3 eVel, double normA, double normB) {
  double waveImp{GetModeImpedance(modeType, n, m, omega)};
  TVector3 j{-TMath::Qe() * eVel};
  ComplexVector3 jComplex{j};

  ComplexVector3 eTrans{GetNormalisedEField(modeType, m, n, ePos, omega)};
  eTrans.SetZ(std::complex<double>{0.0, 0.0});
  ComplexVector3 eAxial{GetNormalisedEField(modeType, m, n, ePos, omega)};
  eAxial.SetX(std::complex<double>{0.0, 0.0});
  eAxial.SetY(std::complex<double>{0.0, 0.0});
  ComplexVector3 subVec{eTrans - eAxial};

  std::complex<double> APlus = (subVec).Dot(jComplex) * (-waveImp / 2.0);
  return APlus;
}

void rad::RectangularWaveguide::CalculatePn(Mode_t modeType, unsigned int n,
                                            unsigned int m, double omega,
                                            unsigned int nSurfPnts) {
  double sum{0};
  // Calculate element of area
  const double elArea{(GetLongDimension() / double(nSurfPnts)) *
                      (GetShortDimension() / double(nSurfPnts))};

  for (int ix{0}; ix < nSurfPnts; ix++) {
    double thisx{-GetLongDimension() / 2.0 +
                 GetLongDimension() / (2.0 * double(nSurfPnts)) +
                 GetLongDimension() * double(ix) / double(nSurfPnts)};
    for (int iy{0}; iy < nSurfPnts; iy++) {
      double thisy{-GetShortDimension() / 2.0 +
                   GetShortDimension() / (2.0 * double(nSurfPnts)) +
                   GetShortDimension() * double(iy) / double(nSurfPnts)};

      TVector3 surfacePos1(thisx, thisy, -0.01);
      TVector3 surfacePos2(thisx, thisy, 0.01);
      ComplexVector3 eTrans1{
          GetNormalisedEField(modeType, n, m, surfacePos1, omega)};
      eTrans1.SetZ(0);
      ComplexVector3 hTrans1{
          GetNormalisedHField(modeType, n, m, surfacePos1, omega)};
      hTrans1.SetZ(0);
      ComplexVector3 eTrans2{
          GetNormalisedEField(modeType, n, m, surfacePos2, omega)};
      eTrans2.SetZ(0);
      ComplexVector3 hTrans2{
          GetNormalisedHField(modeType, n, m, surfacePos2, omega)};
      hTrans2.SetZ(0);
      sum += (eTrans1.Cross(hTrans1)).Z().real() * elArea;
      sum += (eTrans2.Cross(hTrans2)).Z().real() * elArea;
    }
  }
  SetPn(sum * 2);
}