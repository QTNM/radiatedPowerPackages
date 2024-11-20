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

double rad::RectangularWaveguide::GetCutoffWavenumber(WaveguideMode mode) {
  unsigned int m{mode.GetModeIndex1()};
  unsigned int n{mode.GetModeIndex2()};
  // Invalid modes
  if ((m == 0 && n == 0) || m < 0 || n < 0) {
    std::cout << "Invalid rectangular waveguide mode: " << m << ", " << n
              << std::endl;
    return -1;
  } else {
    double k_c{sqrt(pow(double(m) * TMath::Pi() / a, 2) +
                    pow(double(n) * TMath::Pi() / b, 2))};
    return k_c;
  }
}

rad::ComplexVector3 rad::RectangularWaveguide::GetModeEFieldComplex(
    WaveguideMode mode, TVector3 pos, double omega, double A, double B) {
  if (pos.X() > a / 2 || pos.Y() > b / 2 || pos.Z() > d / 2)
    ComplexVector3(0, 0, 0);

  double x{pos.X() + a / 2.0};
  double y{pos.Y() + b / 2.0};
  double z{pos.Z() + d / 2.0};
  std::complex<double> i{0.0, 1.0};
  double k_c{GetCutoffWavenumber(mode)};
  std::complex<double> betaSq{pow(omega / TMath::C(), 2) - k_c * k_c};
  std::complex<double> beta{sqrt(betaSq)};

  ModeType modeType{mode.GetModeType()};
  unsigned int m{mode.GetModeIndex1()};
  unsigned int n{mode.GetModeIndex2()};

  if (modeType == ModeType::kTE) {
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
  } else if (modeType == ModeType::kTM) {
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

TVector3 rad::RectangularWaveguide::GetModeEField(TVector3 pos,
                                                  WaveguideMode mode, double A,
                                                  double omega, bool state) {
  if (pos.X() > a / 2 || pos.Y() > b / 2 || pos.Z() > d / 2)
    return TVector3(0, 0, 0);

  double x{pos.X() + a / 2.0};
  double y{pos.Y() + b / 2.0};
  double z{pos.Z() + d / 2.0};

  double k_c{GetCutoffWavenumber(mode)};
  double m{double(mode.GetModeIndex1())};
  double n{double(mode.GetModeIndex2())};
  double betaSq{pow(omega / TMath::C(), 2) - k_c * k_c};
  double beta{sqrt(betaSq)};

  if (mode.GetModeType() == ModeType::kTE) {
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
  } else if (mode.GetModeType() == ModeType::kTM) {
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

TVector3 rad::RectangularWaveguide::GetNormalisedEField(WaveguideMode mode,
                                                        TVector3 pos,
                                                        double omega) {
  if (pos.X() > a / 2 || pos.Y() > b / 2 || pos.Z() > d / 2)
    return TVector3(0, 0, 0);

  double x{pos.X() + a / 2.0};
  double y{pos.Y() + b / 2.0};
  double z{pos.Z() + d / 2.0};
  double k_c{GetCutoffWavenumber(mode)};
  double betaSq{pow(omega / TMath::C(), 2) - k_c * k_c};
  double beta{sqrt(betaSq)};

  ModeType modeType{mode.GetModeType()};
  unsigned int m{mode.GetModeIndex1()};
  unsigned int n{mode.GetModeIndex2()};

  if (modeType == ModeType::kTE) {
    double Ex{(double(n) / b) * cos(double(m) * TMath::Pi() * x / a) *
              sin(double(n) * TMath::Pi() * y / b)};
    double Ey{-(double(m) / a) * sin(double(m) * TMath::Pi() * x / a) *
              cos(double(n) * TMath::Pi() * y / b)};
    double Ez{0};
    TVector3 eField{Ex, Ey, Ez};
    eField *= sqrt(4 * a * b / (pow(a * double(n), 2) + pow(b * double(m), 2)));
    if (m == 0 || n == 0) eField *= (1.0 / sqrt(2.0));
    return eField;
  } else if (modeType == ModeType::kTM) {
    double Ex{-(double(m) / a) * cos(double(m) * TMath::Pi() * x / a) *
              sin(double(n) * TMath::Pi() * y / b)};
    double Ey{-(double(n) / b) * sin(double(m) * TMath::Pi() * x / a) *
              cos(double(n) * TMath::Pi() * y / b)};
    double Ez{k_c * k_c / M_PI * sin(double(m) * TMath::Pi() * x / a) *
              sin(double(n) * TMath::Pi() * y / b)};
    TVector3 eField{Ex, Ey, Ez};
    eField *= sqrt(4 * a * b / (pow(a * double(m), 2) + pow(b * double(n), 2)));
    return eField;
  } else {
    // We have a TEM mode
    std::cout << "TEM modes are not supported by rectangular waveguides."
              << std::endl;
    TVector3 eField(0, 0, 0);
    return eField;
  }
}

rad::ComplexVector3 rad::RectangularWaveguide::GetModeHFieldComplex(
    WaveguideMode mode, TVector3 pos, double omega, double A, double B) {
  if (pos.X() > a / 2 || pos.Y() > b / 2 || pos.Z() > d / 2)
    ComplexVector3(0, 0, 0);

  double x{pos.X() + a / 2.0};
  double y{pos.Y() + b / 2.0};
  double z{pos.Z() + d / 2.0};
  std::complex<double> i{0.0, 1.0};
  double k_c{GetCutoffWavenumber(mode)};
  std::complex<double> betaSq{pow(omega / TMath::C(), 2) - k_c * k_c};
  std::complex<double> beta{sqrt(betaSq)};

  ModeType modeType{mode.GetModeType()};
  unsigned int m{mode.GetModeIndex1()};
  unsigned int n{mode.GetModeIndex2()};

  if (modeType == ModeType::kTE) {
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
  } else if (modeType == ModeType::kTM) {
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

TVector3 rad::RectangularWaveguide::GetModeHField(WaveguideMode mode,
                                                  TVector3 pos, double omega,
                                                  double A, double B) {
  ComplexVector3 field{GetModeHFieldComplex(mode, pos, omega, A, B)};
  return field.Real();
}

TVector3 rad::RectangularWaveguide::GetModalHField(WaveguideMode mode,
                                                   TVector3 pos, double omega,
                                                   double A) {
  if (pos.X() > a || pos.Y() > b || pos.Z() > d) return TVector3(0, 0, 0);

  double x{pos.X() + a / 2.0};
  double y{pos.Y() + b / 2.0};
  double z{pos.Z() + d / 2.0};
  double k_c{GetCutoffWavenumber(mode)};
  double betaSq{pow(omega / TMath::C(), 2) - k_c * k_c};
  double beta{sqrt(betaSq)};

  ModeType modeType{mode.GetModeType()};
  unsigned int m{mode.GetModeIndex1()};
  unsigned int n{mode.GetModeIndex2()};

  if (modeType == ModeType::kTE) {
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
  } else if (modeType == ModeType::kTM) {
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
    WaveguideMode mode, TVector3 pos, double omega) {
  if (pos.X() > a / 2 || pos.Y() > b / 2 || pos.Z() > d / 2)
    ComplexVector3(0, 0, 0);

  double x{pos.X() + a / 2.0};
  double y{pos.Y() + b / 2.0};
  double z{pos.Z() + d / 2.0};

  double k_c{GetCutoffWavenumber(mode)};
  double betaSq{pow(omega / TMath::C(), 2) - k_c * k_c};
  double beta{sqrt(betaSq)};
  std::complex<double> i{0, 1};
  double waveImp{GetModeImpedance(mode, omega)};
  TVector3 unitZ{0, 0, 1};
  ComplexVector3 unitZComplex{unitZ};

  ModeType modeType{mode.GetModeType()};
  unsigned int m{mode.GetModeIndex1()};
  unsigned int n{mode.GetModeIndex2()};

  if (modeType == ModeType::kTE) {
    ComplexVector3 eFieldTrans{GetNormalisedEField(mode, pos, omega)};
    eFieldTrans.SetZ(std::complex<double>{0, 0});
    ComplexVector3 hFieldTrans{unitZComplex.Cross(eFieldTrans) *
                               (1.0 / waveImp)};
    std::complex<double> Hz{(-2.0 * i * k_c) / (beta * waveImp * sqrt(a * b)) *
                            cos(double(m) * TMath::Pi() * x / a) *
                            cos(double(n) * TMath::Pi() * y / b)};
    ComplexVector3 hField{hFieldTrans.X(), hFieldTrans.Y(), Hz};
    return hField;
  } else if (modeType == ModeType::kTM) {
    ComplexVector3 eFieldTrans{GetNormalisedEField(mode, pos, omega)};
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

double rad::RectangularWaveguide::GetCutoffFrequency(WaveguideMode mode) {
  unsigned int m{mode.GetModeIndex1()};
  unsigned int n{mode.GetModeIndex2()};
  if ((m == 0 && n == 0) || m < 0 || n < 0) {
    std::cout << "Invalid rectangular waveguide mode: " << m << ", " << n
              << std::endl;
    return -1;
  } else {
    double k_c{GetCutoffWavenumber(mode)};
    return k_c * TMath::C() / (2 * TMath::Pi());
  }
}

double rad::RectangularWaveguide::GetEFieldIntegral(WaveguideMode mode,
                                                    double omega, double A,
                                                    int nSurfPnts, bool state) {
  double integral{0};
  const double dX{a / double(nSurfPnts)};
  const double dY{b / double(nSurfPnts)};
  const double area{dX * dY};
  // Loop through points in the field to calculate the integral
  for (unsigned int iX{0}; iX < nSurfPnts; iX++) {
    const double x{(-a / 2) + dX / 2 + double(iX) * dX};
    for (unsigned int iY{0}; iY < nSurfPnts; iY++) {
      const double y{(-b / 2) + dY / 2 + double(iY) * dY};
      TVector3 surfacePos{x, y, 0};
      TVector3 eTrans{GetModeEField(surfacePos, mode, A, omega, true)};
      eTrans.SetZ(0);
      integral += eTrans.Dot(eTrans) * area;
    }
  }

  if (integral == 0)
    return 0;
  else
    return integral;
}

double rad::RectangularWaveguide::GetFieldAmp(WaveguideMode mode, double omega,
                                              TVector3 ePos, TVector3 eVel,
                                              double normA, bool state,
                                              bool isPositive) {
  double waveImp{GetModeImpedance(mode, omega)};
  TVector3 j{-TMath::Qe() * eVel};
  TVector3 jComplex{j};

  TVector3 eTrans{GetNormalisedEField(mode, ePos, omega)};
  eTrans.SetZ(0);
  TVector3 eAxial{GetNormalisedEField(mode, ePos, omega)};
  eAxial.SetX(0);
  eAxial.SetY(0);

  TVector3 eField(0, 0, 0);
  if (isPositive) {
    eField = eTrans - eAxial;
  } else {
    eField = eTrans + eAxial;
  }

  double A{eField.Dot(jComplex) * (-waveImp / 2.0)};
  return A;
}

void rad::RectangularWaveguide::CalculatePn(WaveguideMode mode, double omega,
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
      ComplexVector3 eTrans1{GetNormalisedEField(mode, surfacePos1, omega)};
      eTrans1.SetZ(0);
      ComplexVector3 hTrans1{GetNormalisedHField(mode, surfacePos1, omega)};
      hTrans1.SetZ(0);
      ComplexVector3 eTrans2{GetNormalisedEField(mode, surfacePos2, omega)};
      eTrans2.SetZ(0);
      ComplexVector3 hTrans2{GetNormalisedHField(mode, surfacePos2, omega)};
      hTrans2.SetZ(0);
      sum += (eTrans1.Cross(hTrans1)).Z().real() * elArea;
      sum += (eTrans2.Cross(hTrans2)).Z().real() * elArea;
    }
  }
  SetPn(sum * 2);
}