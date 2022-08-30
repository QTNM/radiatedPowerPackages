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

TVector3 rad::RectangularWaveguide::GetModeEField(Mode_t modeType, int m, int n, TVector3 pos, double omega, double A)
{
  double x{ pos.X() };
  double y{ pos.Y() };
  double z{ pos.Z() };
  std::complex<double> i{ 0.0, 1.0 };
  double k_c{ GetCutoffWavenumber(m, n) };
  std::complex<double> beta{ sqrt(pow(omega/TMath::C(), 2) - k_c*k_c) };
  
  if (modeType == kTE) {
    // We have a TE mode
    std::complex<double> Ex{ (i*omega*MU0*double(n)*TMath::Pi()/(k_c*k_c*b)) * A * cos(double(m)*TMath::Pi()*x/a) * sin(double(n)*TMath::Pi()*y/b) * exp(-i*beta*z) };
    std::complex<double> Ey{ (-1.0*i*omega*MU0*double(m)*TMath::Pi()/(k_c*k_c*a)) * A * sin(double(m)*TMath::Pi()*x/a) * cos(double(n)*TMath::Pi()*y/b) * exp(-i*beta*z) };
    std::complex<double> Ez{ 0.0, 0.0 };
    TVector3 eField{ Ex.real(), Ey.real(), Ez.real() };
    return eField;
  }
  else if (modeType == kTM) {
    // We have a TM mode
    std::complex<double> Ex{ (-1.0*i*beta*double(m)*TMath::Pi()/(k_c*k_c*a)) * A * cos(double(m)*TMath::Pi()*x/a) * sin(double(n)*TMath::Pi()*y/b) * exp(-i*beta*z) };
    std::complex<double> Ey{ (-1.0*i*beta*double(n)*TMath::Pi()/(k_c*k_c*b)) * A * sin(double(m)*TMath::Pi()*x/a) * cos(double(n)*TMath::Pi()*y/b) * exp(-i*beta*z) };
    std::complex<double> Ez{ A * sin(double(m)*TMath::Pi()*x/a) * sin(double(n)*TMath::Pi()*y/b) * exp(-i*beta*z) };
    TVector3 eField{ Ex.real(), Ey.real(), Ez.real() };
    return eField;
  }
  else {
    // We have a TEM mode
    std::cout<<"TEM modes are not supported by circular waveguides."<<std::endl;
    TVector3 eField{ 0, 0, 0 };
    return eField;
  }
}

TVector3 rad::RectangularWaveguide::GetModeHField(Mode_t modeType, int m, int n, TVector3 pos, double omega, double A)
{
  double x{ pos.X() };
  double y{ pos.Y() };
  double z{ pos.Z() };
  std::complex<double> i{ 0.0, 1.0 };
  double k_c{ GetCutoffWavenumber(m, n) };
  std::complex<double> beta{ sqrt(pow(omega/TMath::C(), 2) - k_c*k_c) };

  if (modeType == kTE) {
    // We have a TE mode
    std::complex<double> Hx{ (i*beta*double(m)*TMath::Pi()/(k_c*k_c*a)) * A * sin(double(m)*TMath::Pi()*x/a) * cos(double(n)*TMath::Pi()*y/b) * exp(-i*beta*z) };
    std::complex<double> Hy{ (i*beta*double(n)*TMath::Pi()/(k_c*k_c*b)) * A * cos(double(m)*TMath::Pi()*x/a) * sin(double(n)*TMath::Pi()*y/b) * exp(-i*beta*z) };
    std::complex<double> Hz{ A * cos(double(m)*TMath::Pi()*x/a) * cos(double(n)*TMath::Pi()*y/b) * exp(-i*beta*z) };
    TVector3 hField{ Hx.real(), Hy.real(), Hz.real() };
    return hField;
  }
  else if (modeType == kTM) {
    // We have a TM mode
    std::complex<double> Hx{ (i*omega*EPSILON0*double(n)*TMath::Pi()/(k_c*k_c*b)) * A * sin(double(m)*TMath::Pi()*x/a) * cos(double(n)*TMath::Pi()*y/b) * exp(-i*beta*z) };
    std::complex<double> Hy{ (-1.0*i*omega*EPSILON0*double(m)*TMath::Pi()/(k_c*k_c*a)) * A * cos(double(m)*TMath::Pi()*x/a) * sin(double(n)*TMath::Pi()*y/b) * exp(-i*beta*z) };
    std::complex<double> Hz{ 0.0, 0.0 };
    TVector3 hField{ Hx.real(), Hy.real(), Hz.real() };
    return hField;
  }
  else {
    // We have a TEM mode
    std::cout<<"TEM modes are not supported by circular waveguides."<<std::endl;
    TVector3 hField{ 0, 0, 0 };
    return hField;
  }
}
