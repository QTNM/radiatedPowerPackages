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

TVector3 rad::RectangularWaveguide::GetModeEField(TVector3 pos,
                                                  WaveguideMode mode, double A,
                                                  double omega, bool state) {
  if (omega < GetCutoffFrequency(mode)) {
    std::cout
        << "Mode is below cutoff frequency and will not propagate. Returning 0."
        << std::endl;
    return TVector3(0, 0, 0);
  }
  double x{pos.X() + a / 2.0};
  double y{pos.Y() + b / 2.0};
  double z{pos.Z() + d / 2.0};
  if (x > a || x < 0 || y > b || y < 0) return TVector3(0, 0, 0);

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

TVector3 rad::RectangularWaveguide::GetModeHField(TVector3 pos,
                                                  WaveguideMode mode, double A,
                                                  double omega, bool state) {
  if (omega < GetCutoffFrequency(mode)) {
    std::cout
        << "Mode is below cutoff frequency and will not propagate. Returning 0."
        << std::endl;
    return TVector3(0, 0, 0);
  }
  double x{pos.X() + a / 2.0};
  double y{pos.Y() + b / 2.0};
  double z{pos.Z() + d / 2.0};
  if (x > a || x < 0 || y > b || y < 0) return TVector3(0, 0, 0);

  double k_c{GetCutoffWavenumber(mode)};
  double m{double(mode.GetModeIndex1())};
  double n{double(mode.GetModeIndex2())};
  double betaSq{pow(omega / TMath::C(), 2) - k_c * k_c};
  double beta{sqrt(betaSq)};

  double HX{0};
  double HY{0};
  double HZ{0};
  if (mode.GetModeType() == ModeType::kTE) {
    HX = A * beta * double(m) * M_PI * sin(double(m) * M_PI * x / a) *
         cos(double(n) * M_PI * y / b) / (k_c * k_c * a);
    HY = A * beta * double(n) * M_PI * cos(double(m) * M_PI * x / a) *
         sin(double(n) * M_PI * y / b) / (k_c * k_c * b);
    HZ = A * cos(double(m) * M_PI * x / a) * cos(double(n) * M_PI * y / b);
  } else if (mode.GetModeType() == ModeType::kTM) {
    HX = A * omega * EPSILON0 * double(n) * M_PI *
         sin(double(m) * M_PI * x / a) * cos(double(n) * M_PI * y / b) /
         (k_c * k_c * b);
    HY = -A * omega * EPSILON0 * double(m) * M_PI *
         cos(double(m) * M_PI * x / a) * sin(double(n) * M_PI * y / b) /
         (k_c * k_c * a);
  } else {
    std::cout << "TEM modes are not supported by rectangular waveguides."
              << std::endl;
    return TVector3(0, 0, 0);
  }
  return TVector3(HX, HY, HZ);
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

      TVector3 surfacePos(thisx, thisy, 0);
      TVector3 eTrans{GetModeEField(surfacePos, mode, 1, omega, true)};
      eTrans.SetZ(0);
      TVector3 hTrans{GetModeHField(surfacePos, mode, 1, omega, true)};
      hTrans.SetZ(0);
      sum += (eTrans.Cross(hTrans)).Dot(TVector3(0, 0, 1)) * elArea;
    }
  }
  SetPn(sum * 2);
}