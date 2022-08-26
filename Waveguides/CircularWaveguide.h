/*
  CircularWaveguide.h

  Class describing a cylindrical waveguide and the modes within
*/

#ifndef CIRCULAR_WAVEGUIDE_H
#define CIRCULAR_WAVEGUIDE_H

#include "TVector3.h"

namespace rad {

  class CircularWaveguide {

  private:
    double a; // Waveguide inner radius (in metres)
    double d; // Waveguide length (in metres)
    // double r_s; // Surface resistance
    
  public:
    enum Mode_t {
      kTE, kTM
    };
    
    /// Parametrised  constructor
    /// \param radius Inner radius of the waveguide (in metres)
    /// \param length Length of the waveguide (in metres)
    CircularWaveguide(double radius, double length) : a{radius}, d{radius} {}

    /// Return the inner radius of the waveguide
    /// \Returns The waveguide inner radius (in metres)
    double GetInnerRadius() { return a; }

    /// Return the length of the waveguide
    /// \Returns The waveguide length (in metres)
    double GetLength() { return d; }

    /// Gets the electric field vector for a given mode at a point
    /// \param modeType The mode type to get (either TE or TM)
    /// \param n The angular number of the mode
    /// \param m The radial number of the mode
    /// \param rho Radial position in cylindrical coordinates in metres
    /// \param phi Angular position in cylindrical coordinates in radians
    /// \param z Longitudinal position in metres
    /// \param omega Angular frequency of the chosen wave
    /// \param A Arbitrary amplitude for part of solution (default = 1)
    /// \param B Arbitrary amplitude for part of solution (default = 0)
    /// \Returns The mode electric field vector at the supplied point
    TVector3 GetModeEField(Mode_t modeType, int n, int m, double rho, double phi, double z, double omega, double A=1, double B=0);

    /// Gets the electric field vector for a given mode at a point
    /// \param modeType The mode type to get (either TE or TM)
    /// \param n The angular number of the mode
    /// \param m The radial number of the mode
    /// \param pos The position vector (in metres)
    /// \param omega Angular frequency of the chosen wave
    /// \param A Arbitrary amplitude for part of solution (default = 1)
    /// \param B Arbitrary amplitude for part of solution (default = 0)
    /// \Returns The mode electric field vector at the supplied point
    TVector3 GetModeEField(Mode_t modeType, int n, int m, TVector3 pos, double omega, double A=1, double B=0);

    /// Gets the H field vector for a given mode at a point
    /// \param modeType The mode type to get (either TE or TM)
    /// \param n The angular number of the mode
    /// \param m The radial number of the mode
    /// \param rho Radial position in cylindrical coordinates in metres
    /// \param phi Angular position in cylindrical coordinates in radians
    /// \param z Longitudinal position in metres
    /// \param omega Angular frequency of the chosen wave
    /// \param A Arbitrary amplitude for part of solution (default = 1)
    /// \param B Arbitrary amplitude for part of solution (default = 0)
    /// \Returns The mode H field vector at the supplied point
    TVector3 GetModeHField(Mode_t modeType, int n, int m, double rho, double phi, double z, double omega, double A=1, double B=0);

    /// Gets the H field vector for a given mode at a point
    /// \param modeType The mode type to get (either TE or TM)
    /// \param n The angular number of the mode
    /// \param m The radial number of the mode
    /// \param pos The position vector (in metres)
    /// \param omega Angular frequency of the chosen wave
    /// \param A Arbitrary amplitude for part of solution (default = 1)
    /// \param B Arbitrary amplitude for part of solution (default = 0)
    /// \Returns The mode H field vector at the supplied point
    TVector3 GetModeHField(Mode_t modeType, int n, int m, TVector3 pos, double omega, double A=1, double B=0);
    
  };
}

#endif
