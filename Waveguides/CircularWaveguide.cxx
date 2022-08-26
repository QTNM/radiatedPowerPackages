/// CircularWaveguide.cxx

#include "Waveguides/CircularWaveguide.h"
#include "BasicFunctions/BasicFunctions.h"

#include <boost/math/special_functions/bessel.hpp>
#include <boost/math/special_functions/bessel_prime.hpp>

#include <complex>
#include <cmath>

#include "TVector3.h"
#include "TMath.h"

TVector3 rad::CircularWaveguide::GetModeEField(Mode_t modeType, int n, int m, double rho, double phi, double z, double omega, double A, double B)
{
  std::complex<double> i{ 0.0, 1.0 };
  if (modeType == kTE) {
    // We have a TE mode
    double pnmPrime{ GetBesselPrimeZero(n, m) }; 
    double k_c{ pnmPrime / a };
    std::complex<double> beta{ sqrt(pow(omega/TMath::C(), 2) - k_c*k_c) };

    std::complex<double> Ez{ 0.0, 0.0 };
    std::complex<double> Erho{ (-1.0*omega*MU0*double(n)/(k_c*k_c*rho)) * ( A*cos(double(n)*phi) - B*sin(double(n)*phi) ) * boost::math::cyl_bessel_j(n, k_c*rho) };
    Erho *= i * exp( -1.0*i*beta*z );
    std::complex<double> Ephi{ (omega*MU0/k_c)*( A*sin(double(n)*phi) + B*cos(double(n)*phi) ) * boost::math::cyl_bessel_j_prime(n, k_c*rho) };
    Ephi *= i * exp( -1.0*i*beta*z );
    
    TVector3 eField{ Erho.real()*cos(phi) - Ephi.real()*sin(phi), Erho.real()*sin(phi) + Ephi.real()*cos(phi), Ez.real() };
    return eField;
  }
  else {
    // We have a TM mode
    double pnm{ boost::math::cyl_bessel_j_zero(double(n), m) };
    double k_c{ pnm / a };
    std::complex<double> beta{ sqrt(pow(omega/TMath::C(), 2) - k_c*k_c) };

    std::complex<double> Ez{ (A*sin(double(n)*phi) + B*cos(double(n)*phi)) * boost::math::cyl_bessel_j(n, k_c*rho), 0.0 };
    Ez *= exp( -1.0*i*beta*z );
    std::complex<double> Erho{ (-1.0*beta/k_c) * (A*sin(double(n)*phi) + B*cos(double(n)*phi)) * boost::math::cyl_bessel_j_prime(n, k_c*rho) };
    Erho *= i * exp( -1.0*i*beta*z );
    std::complex<double> Ephi{ (-1.0*beta*double(n)/(k_c*k_c*rho)) * (A*cos(double(n)*phi) - B*sin(double(n)*phi)) * boost::math::cyl_bessel_j(n, k_c*rho) };
    Ephi *= i * exp( -1.0*i*beta*z );

    TVector3 eField{ Erho.real()*cos(phi) - Ephi.real()*sin(phi), Erho.real()*sin(phi) + Ephi.real()*cos(phi), Ez.real() };
    return eField;    
  }
}

TVector3 rad::CircularWaveguide::GetModeEField(Mode_t modeType, int n, int m, TVector3 pos, double omega, double A, double B)
{
  return GetModeEField(modeType, n, m, pos.Perp(), pos.Phi(), pos.Z(), omega, A, B);
}

TVector3 rad::CircularWaveguide::GetModeHField(Mode_t modeType, int n, int m, double rho, double phi, double z, double omega, double A, double B)
{
  std::complex<double> i{ 0.0, 1.0 };
  if (modeType == kTE) {
    // We have a TE mode
    double pnmPrime{ GetBesselPrimeZero(n, m) }; 
    double k_c{ pnmPrime / a };
    std::complex<double> beta{ sqrt(pow(omega/TMath::C(), 2) - k_c*k_c) };

    std::complex<double> Hz{ (A*sin(double(n)*phi) + B*cos(double(n)*phi)) * boost::math::cyl_bessel_j(n, k_c*rho), 0.0 };
    Hz *= exp( -1.0*i*beta*z );
    std::complex<double> Hrho{ (-1.0*beta/k_c) * (A*sin(double(n)*phi) + B*cos(double(n)*phi)) * boost::math::cyl_bessel_j_prime(n, k_c*rho) };
    Hrho *= i * exp( -1.0*i*beta*z );
    std::complex<double> Hphi{ (-1.0*beta*double(n)/(k_c*k_c*rho)) * (A*cos(double(n)*phi) - B*sin(double(n)*phi)) * boost::math::cyl_bessel_j(n, k_c*rho) };
    Hphi *= i * exp( -1.0*i*beta*z );

    TVector3 hField{ Hrho.real()*cos(phi) - Hphi.real()*sin(phi), Hrho.real()*sin(phi) + Hphi.real()*cos(phi), Hz.real() };
    return hField; 
  }
  else {
    // We have a TM mode
    double pnm{ boost::math::cyl_bessel_j_zero(double(n), m) };
    double k_c{ pnm / a };
    std::complex<double> beta{ sqrt(pow(omega/TMath::C(), 2) - k_c*k_c) };

    std::complex<double> Hz{ 0.0, 0.0 };
    std::complex<double> Hrho{ (omega*EPSILON0*double(n)/(k_c*k_c*rho)) * (A*cos(double(n)*phi) - B*sin(double(n)*phi)) * boost::math::cyl_bessel_j(n, k_c*rho) };
    Hrho *= i * exp( -1.0*i*beta*z );
    std::complex<double> Hphi{ (-1.0*omega*EPSILON0/k_c) * (A*sin(double(n)*phi) + B*cos(double(n)*phi)) * boost::math::cyl_bessel_j_prime(n, k_c*rho) };
    Hphi *= i * exp( -1.0*i*beta*z );

    TVector3 hField{ Hrho.real()*cos(phi) - Hphi.real()*sin(phi), Hrho.real()*sin(phi) + Hphi.real()*cos(phi), Hz.real() };
    return hField; 
  }
}

TVector3 rad::CircularWaveguide::GetModeHField(Mode_t modeType, int n, int m, TVector3 pos, double omega, double A, double B)
{
  return GetModeHField(modeType, n, m, pos.Perp(), pos.Phi(), pos.Z(), omega, A, B);
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
  else {
    double pnm{ boost::math::cyl_bessel_j_zero(double(n), m) };
    double k_c{ pnm / a };
    double k{ omega/TMath::C() };
    double beta{ sqrt(pow(k*k, 2) - k_c*k_c) };
    return beta * sqrt(MU0/EPSILON0) / k;
  }
}
