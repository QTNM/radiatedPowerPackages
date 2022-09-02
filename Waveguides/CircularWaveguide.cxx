/// CircularWaveguide.cxx

#include "Waveguides/CircularWaveguide.h"
#include "BasicFunctions/BasicFunctions.h"

#include <boost/math/special_functions/bessel.hpp>
#include <boost/math/special_functions/bessel_prime.hpp>

#include <complex>
#include <cmath>
#include <iostream>

#include "TVector3.h"
#include "TMath.h"

rad::ComplexVector3 rad::CircularWaveguide::GetModeEFieldComplex(Mode_t modeType, int n, int m, TVector3 pos, double omega, double A, double B)
{
  double rho{ pos.Perp() };
  double phi{ pos.Phi() };
  double z{ pos.Z() + d/2.0 };
  std::complex<double> i{ 0.0, 1.0 };
  
  if (modeType == kTE) {
    // We have a TE mode
    double pnmPrime{ GetBesselPrimeZero(n, m) }; 
    double k_c{ pnmPrime / a };
    std::complex<double> betaSq{ pow(omega/TMath::C(), 2) - k_c*k_c };
    std::complex<double> beta{ sqrt(betaSq) };
    
    std::complex<double> Ez{ 0.0, 0.0 };
    std::complex<double> Erho{ (-1.0*omega*MU0*double(n)/(k_c*k_c*rho)) * ( A*cos(double(n)*phi) - B*sin(double(n)*phi) ) * boost::math::cyl_bessel_j(n, k_c*rho) };
    Erho *= i * exp( -1.0*i*beta*z );
    std::complex<double> Ephi{ (omega*MU0/k_c)*( A*sin(double(n)*phi) + B*cos(double(n)*phi) ) * boost::math::cyl_bessel_j_prime(n, k_c*rho) };
    Ephi *= i * exp( -1.0*i*beta*z );
    
    ComplexVector3 eField{ Erho*cos(phi) - Ephi*sin(phi), Erho*sin(phi) + Ephi*cos(phi), Ez };
    return eField;
  }
  else if (modeType == kTM) {
    // We have a TM mode
    double pnm{ boost::math::cyl_bessel_j_zero(double(n), m) };
    double k_c{ pnm / a };
    std::complex<double> betaSq{ pow(omega/TMath::C(), 2) - k_c*k_c };
    std::complex<double> beta{ sqrt(betaSq) };
      
    std::complex<double> Ez{ (A*sin(double(n)*phi) + B*cos(double(n)*phi)) * boost::math::cyl_bessel_j(n, k_c*rho), 0.0 };
    Ez *= exp( -1.0*i*beta*z );
    std::complex<double> Erho{ (-1.0*beta/k_c) * (A*sin(double(n)*phi) + B*cos(double(n)*phi)) * boost::math::cyl_bessel_j_prime(n, k_c*rho) };
    Erho *= i * exp( -1.0*i*beta*z );
    std::complex<double> Ephi{ (-1.0*beta*double(n)/(k_c*k_c*rho)) * (A*cos(double(n)*phi) - B*sin(double(n)*phi)) * boost::math::cyl_bessel_j(n, k_c*rho) };
    Ephi *= i * exp( -1.0*i*beta*z );

    ComplexVector3 eField{ Erho*cos(phi) - Ephi*sin(phi), Erho*sin(phi) + Ephi*cos(phi), Ez };
    return eField;    
  }
  else {
    // This is a TEM mode that which is not supported by this waveguide
    std::cout<<"TEM modes are not supported by circular waveguides."<<std::endl;
    ComplexVector3 eField{ std::complex<double>{0}, std::complex<double>{0}, std::complex<double>{0} };
    return eField;
  }
}

TVector3 rad::CircularWaveguide::GetModeEField(Mode_t modeType, int n, int m, TVector3 pos, double omega, double A, double B)
{
  ComplexVector3 field{ GetModeEFieldComplex(modeType, n, m, pos, omega, A, B) };
  return field.Real();
}

rad::ComplexVector3 rad::CircularWaveguide::GetModeHFieldComplex(Mode_t modeType, int n, int m, TVector3 pos, double omega, double A, double B)
{
  double rho{ pos.Perp() };
  double phi{ pos.Phi() };
  double z{ pos.Z() };
  std::complex<double> i{ 0.0, 1.0 };
  
  if (modeType == kTE) {
    // We have a TE mode
    double pnmPrime{ GetBesselPrimeZero(n, m) }; 
    double k_c{ pnmPrime / a };
    std::complex<double> betaSq{ pow(omega/TMath::C(), 2) - k_c*k_c };
    std::complex<double> beta{ sqrt(betaSq) };
    
    std::complex<double> Hz{ (A*sin(double(n)*phi) + B*cos(double(n)*phi)) * boost::math::cyl_bessel_j(n, k_c*rho), 0.0 };
    Hz *= exp( -1.0*i*beta*z );
    std::complex<double> Hrho{ (-1.0*beta/k_c) * (A*sin(double(n)*phi) + B*cos(double(n)*phi)) * boost::math::cyl_bessel_j_prime(n, k_c*rho) };
    Hrho *= i * exp( -1.0*i*beta*z );
    std::complex<double> Hphi{ (-1.0*beta*double(n)/(k_c*k_c*rho)) * (A*cos(double(n)*phi) - B*sin(double(n)*phi)) * boost::math::cyl_bessel_j(n, k_c*rho) };
    Hphi *= i * exp( -1.0*i*beta*z );

    ComplexVector3 hField{ Hrho*cos(phi) - Hphi*sin(phi), Hrho*sin(phi) + Hphi*cos(phi), Hz };
    return hField; 
  }
  else if (modeType == kTM) {
    // We have a TM mode
    double pnm{ boost::math::cyl_bessel_j_zero(double(n), m) };
    double k_c{ pnm / a };
    std::complex<double> betaSq{ pow(omega/TMath::C(), 2) - k_c*k_c };
    std::complex<double> beta{ sqrt(betaSq) };
    
    std::complex<double> Hz{ 0.0, 0.0 };
    std::complex<double> Hrho{ (omega*EPSILON0*double(n)/(k_c*k_c*rho)) * (A*cos(double(n)*phi) - B*sin(double(n)*phi)) * boost::math::cyl_bessel_j(n, k_c*rho) };
    Hrho *= i * exp( -1.0*i*beta*z );
    std::complex<double> Hphi{ (-1.0*omega*EPSILON0/k_c) * (A*sin(double(n)*phi) + B*cos(double(n)*phi)) * boost::math::cyl_bessel_j_prime(n, k_c*rho) };
    Hphi *= i * exp( -1.0*i*beta*z );

    ComplexVector3 hField{ Hrho*cos(phi) - Hphi*sin(phi), Hrho*sin(phi) + Hphi*cos(phi), Hz };
    return hField; 
  }
  else {
    // This is a TEM mode that which is not supported by this waveguide
    std::cout<<"TEM modes are not supported by circular waveguides."<<std::endl;
    ComplexVector3 eField{ std::complex<double>{0}, std::complex<double>{0}, std::complex<double>{0} };
    return eField;
  }
}

TVector3 rad::CircularWaveguide::GetModeHField(Mode_t modeType, int n, int m, TVector3 pos, double omega, double A, double B)
{
  ComplexVector3 field{ GetModeHFieldComplex(modeType, n, m, pos, omega, A, B) };
  return field.Real();
}

double rad::CircularWaveguide::GetModeImpedance(Mode_t modeType, int n, int m, double omega)
{
  if (modeType == kTE) {
    double pnmPrime{ GetBesselPrimeZero(n, m) }; 
    double k_c{ pnmPrime / a };
    double k{ omega/TMath::C() };
    double beta{ sqrt(pow(k*k, 2) - k_c*k_c) };
    return k * sqrt(MU0/EPSILON0) / beta;
  }
  else if (modeType == kTM) {
    double pnm{ boost::math::cyl_bessel_j_zero(double(n), m) };
    double k_c{ pnm / a };
    double k{ omega/TMath::C() };
    double beta{ sqrt(pow(k*k, 2) - k_c*k_c) };
    return beta * sqrt(MU0/EPSILON0) / k;
  }
  else {
    // This is a TEM mode that which is not supported by this waveguide
    std::cout<<"TEM modes are not supported by circular waveguides."<<std::endl;
    return -1;
  }
}

double rad::CircularWaveguide::GetCutoffFrequency(Mode_t modeType, int n, int m)
{
  if (modeType == kTE) {
    double pnmPrime{ GetBesselPrimeZero(n, m) };
    double k_c{ pnmPrime / a };
    double f_c{ k_c * TMath::C() / (2 * TMath::Pi()) };
    return f_c;
  }
  else if (modeType == kTM) {
    double pnm{ boost::math::cyl_bessel_j_zero(double(n), m) };
    double k_c{ pnm / a };
    double f_c{ k_c * TMath::C() / (2 * TMath::Pi()) };
    return f_c;
  }
  else {
    return -1;
  }
}

double rad::CircularWaveguide::GetResonantModeFrequency(Mode_t modeType, int n, int m, int l)
{
  double premult{ TMath::C()/(2*TMath::C()) };
  if (modeType == kTE) {
    double pnmPrime{ GetBesselPrimeZero(n, m) };
    double freq{ premult*sqrt( pow(pnmPrime/a, 2) + pow(double(l)*TMath::Pi()/d, 2) ) };
    return freq;
  }
  else if (modeType == kTM) {
    double pnm{ boost::math::cyl_bessel_j_zero(double(n), m) };
    double freq{ premult*sqrt( pow(pnm/a, 2) + pow(double(l)*TMath::Pi()/d, 2) ) };
    return freq;
  }
  else {
    return -1;
  }
}
