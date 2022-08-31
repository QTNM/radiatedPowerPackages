/// RectangularWaveguide.cxx

#include "Waveguides/RectangularWaveguide.h"
#include "BasicFunctions/Constants.h"

#include <complex>
#include <cmath>
#include <iostream>

#include "TMath.h"
#include "TVector3.h"

rad::RectangularWaveguide::RectangularWaveguide(double longSide, double shortSide)
{
  if (longSide >= shortSide) {
    a = longSide;
    b = shortSide;
  }
  else {
    a = shortSide;
    b = longSide;
  }
}

double rad::RectangularWaveguide::GetCutoffWavenumber(unsigned int m, unsigned int n)
{
  double k_c{ sqrt( pow(double(m)*TMath::Pi()/a, 2) + pow(double(n)*TMath::Pi()/b, 2) ) };
  return k_c;
}

rad::ComplexVector3 rad::RectangularWaveguide::GetModeEFieldComplex(Mode_t modeType, int m, int n, TVector3 pos, double omega, double A, double B)
{
  double x{ pos.X() };
  double y{ pos.Y() };
  double z{ pos.Z() };
  std::complex<double> i{ 0.0, 1.0 };
  double k_c{ GetCutoffWavenumber(m, n) };
  std::complex<double> betaSq{ pow(omega/TMath::C(), 2) - k_c*k_c };
  std::complex<double> beta{ sqrt(betaSq) };
  
  if (modeType == kTE) {
    // We have a TE mode
    std::complex<double> Ex{ (i*omega*MU0*double(n)*TMath::Pi()/(k_c*k_c*b)) * A * cos(double(m)*TMath::Pi()*x/a) * sin(double(n)*TMath::Pi()*y/b) * exp(-i*beta*z) };
    std::complex<double> Ey{ (-1.0*i*omega*MU0*double(m)*TMath::Pi()/(k_c*k_c*a)) * A * sin(double(m)*TMath::Pi()*x/a) * cos(double(n)*TMath::Pi()*y/b) * exp(-i*beta*z) };
    std::complex<double> Ez{ 0.0, 0.0 };
    ComplexVector3 eField{ Ex, Ey, Ez };
    return eField;
  }
  else if (modeType == kTM) {
    // We have a TM mode
    std::complex<double> Ex{ (-1.0*i*beta*double(m)*TMath::Pi()/(k_c*k_c*a)) * A * cos(double(m)*TMath::Pi()*x/a) * sin(double(n)*TMath::Pi()*y/b) * exp(-i*beta*z) };
    std::complex<double> Ey{ (-1.0*i*beta*double(n)*TMath::Pi()/(k_c*k_c*b)) * A * sin(double(m)*TMath::Pi()*x/a) * cos(double(n)*TMath::Pi()*y/b) * exp(-i*beta*z) };
    std::complex<double> Ez{ A * sin(double(m)*TMath::Pi()*x/a) * sin(double(n)*TMath::Pi()*y/b) * exp(-i*beta*z) };
    ComplexVector3 eField{ Ex, Ey, Ez };
    return eField;
  }
  else {
    // We have a TEM mode
    std::cout<<"TEM modes are not supported by circular waveguides."<<std::endl;
    ComplexVector3 eField{ std::complex<double>{0}, std::complex<double>{0}, std::complex<double>{0} };
    return eField;
  }
}

TVector3 rad::RectangularWaveguide::GetModeEField(Mode_t modeType, int m, int n, TVector3 pos, double omega, double A, double B)
{
  ComplexVector3 field{ GetModeEFieldComplex(modeType, m, n, pos, omega, A, B) };
  return field.Real();
}

rad::ComplexVector3 rad::RectangularWaveguide::GetModeHFieldComplex(Mode_t modeType, int m, int n, TVector3 pos, double omega, double A, double B)
{
  double x{ pos.X() };
  double y{ pos.Y() };
  double z{ pos.Z() };
  std::complex<double> i{ 0.0, 1.0 };
  double k_c{ GetCutoffWavenumber(m, n) };
  std::complex<double> betaSq{ pow(omega/TMath::C(), 2) - k_c*k_c };
  std::complex<double> beta{ sqrt(betaSq) };

  if (modeType == kTE) {
    // We have a TE mode
    std::complex<double> Hx{ (i*beta*double(m)*TMath::Pi()/(k_c*k_c*a)) * A * sin(double(m)*TMath::Pi()*x/a) * cos(double(n)*TMath::Pi()*y/b) * exp(-i*beta*z) };
    std::complex<double> Hy{ (i*beta*double(n)*TMath::Pi()/(k_c*k_c*b)) * A * cos(double(m)*TMath::Pi()*x/a) * sin(double(n)*TMath::Pi()*y/b) * exp(-i*beta*z) };
    std::complex<double> Hz{ A * cos(double(m)*TMath::Pi()*x/a) * cos(double(n)*TMath::Pi()*y/b) * exp(-i*beta*z) };
    ComplexVector3 hField{ Hx, Hy, Hz };
    return hField;
  }
  else if (modeType == kTM) {
    // We have a TM mode
    std::complex<double> Hx{ (i*omega*EPSILON0*double(n)*TMath::Pi()/(k_c*k_c*b)) * A * sin(double(m)*TMath::Pi()*x/a) * cos(double(n)*TMath::Pi()*y/b) * exp(-i*beta*z) };
    std::complex<double> Hy{ (-1.0*i*omega*EPSILON0*double(m)*TMath::Pi()/(k_c*k_c*a)) * A * cos(double(m)*TMath::Pi()*x/a) * sin(double(n)*TMath::Pi()*y/b) * exp(-i*beta*z) };
    std::complex<double> Hz{ 0.0, 0.0 };
    ComplexVector3 hField{ Hx, Hy, Hz };
    return hField;
  }
  else {
    // We have a TEM mode
    std::cout<<"TEM modes are not supported by circular waveguides."<<std::endl;
    ComplexVector3 hField{ std::complex<double>{0}, std::complex<double>{0}, std::complex<double>{0} };
    return hField;
  }
}

TVector3 rad::RectangularWaveguide::GetModeHField(Mode_t modeType, int m, int n, TVector3 pos, double omega, double A, double B)
{
  ComplexVector3 field{ GetModeHFieldComplex(modeType, m, n, pos, omega, A, B) };
  return field.Real();
}

double rad::RectangularWaveguide::GetModeImpedance(Mode_t modeType, int m, int n, double omega)
{
  double k_c{ GetCutoffWavenumber(m, n) };
  double k{ omega / TMath::C() };
  double beta{ sqrt(k*k - k_c*k_c) };

  if (modeType == kTE) {
    return k * sqrt(MU0/EPSILON0) / beta;
  }
  else if (modeType == kTM) {
    return beta * sqrt(MU0/EPSILON0) / k;
  }
  else {
    // We have a TEM mode
    std::cout<<"TEM modes are not supported by circular waveguides."<<std::endl;
    return 0;
  }
}

double rad::RectangularWaveguide::GetCutoffFrequency(Mode_t modeType, int m, int n)
{
  double k_c{ GetCutoffWavenumber(m, n) };
  return k_c*TMath::C() / (2*TMath::Pi());
}
