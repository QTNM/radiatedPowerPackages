/*
  CircularWaveguide.h

  Class describing a cylindrical waveguide and the modes within
*/

#ifndef CIRCULAR_WAVEGUIDE_H
#define CIRCULAR_WAVEGUIDE_H

#include "TVector3.h"
#include "Waveguides/IWaveguide.h"

namespace rad {

/// Class describing a circular cylindrical waveguide
/// The centre of the waveguide (in the axial direction) is located at z=0
/// Centre of the circular cross-section located at x = y = 0;
class CircularWaveguide : public IWaveguide {
 private:
  double a;  // Waveguide inner radius (in metres)
  double d;  // Waveguide length (in metres)
             // double r_s; // Surface resistance

 public:
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

  /// Complex electric field vector
  /// \param modeType The mode type to get (either TE or TM)
  /// \param n The angular number of the mode
  /// \param m The radial number of the mode
  /// \param pos The position vector (in metres)
  /// \param omega Angular frequency of the chosen wave
  /// \param A Arbitrary amplitude for part of solution (default = 1)
  /// \param B Arbitrary amplitude for part of solution (default = 0)
  /// \Returns The mode complex electric field vector at the supplied point
  ComplexVector3 GetModeEFieldComplex(Mode_t modeType, int n, int m,
                                      TVector3 pos, double omega, double A = 1,
                                      double B = 0) override;

  /// Gets the real electric field vector for a given mode at a point
  /// \param modeType The mode type to get (either TE or TM)
  /// \param n The angular number of the mode
  /// \param m The radial number of the mode
  /// \param pos The position vector (in metres)
  /// \param omega Angular frequency of the chosen wave
  /// \param A Arbitrary amplitude for part of solution (default = 1)
  /// \param B Arbitrary amplitude for part of solution (default = 0)
  /// \Returns The mode electric field vector at the supplied point
  TVector3 GetModeEField(Mode_t modeType, int n, int m, TVector3 pos,
                         double omega, double A = 1, double B = 0) override;

  ComplexVector3 GetModalEField(Mode_t modeType, int n, int m, TVector3 pos,
                                double omega, double A = 1, double B = 0);

  /// Gets the complex magnetic field strength vector for a given mode at a
  /// point \param modeType The mode type to get (either TE or TM) \param n The
  /// angular number of the mode \param m The radial number of the mode \param
  /// pos The position vector (in metres) \param omega Angular frequency of the
  /// chosen wave \param A Arbitrary amplitude for part of solution (default =
  /// 1) \param B Arbitrary amplitude for part of solution (default = 0)
  /// \Returns The mode H field vector at the supplied point
  ComplexVector3 GetModeHFieldComplex(Mode_t modeType, int n, int m,
                                      TVector3 pos, double omega, double A = 1,
                                      double B = 0) override;

  /// Gets the H field vector for a given mode at a point
  /// \param modeType The mode type to get (either TE or TM)
  /// \param n The angular number of the mode
  /// \param m The radial number of the mode
  /// \param pos The position vector (in metres)
  /// \param omega Angular frequency of the chosen wave
  /// \param A Arbitrary amplitude for part of solution (default = 1)
  /// \param B Arbitrary amplitude for part of solution (default = 0)
  /// \Returns The mode H field vector at the supplied point
  TVector3 GetModeHField(Mode_t modeType, int n, int m, TVector3 pos,
                         double omega, double A = 1, double B = 0) override;

  ComplexVector3 GetModalHField(Mode_t modeType, int n, int m, TVector3 pos,
                                double omega, double A = 1, double B = 0);

  /// Gets the cutoff frequency for a particular waveguide mode
  /// \param modeType The mode type to get (TE, TM, TEM)
  /// \param n The angular number of the mode
  /// \param m The radial number of the mode
  /// \Returns The cutoff frequency in Hertz
  double GetCutoffFrequency(Mode_t modeType, int n, int m) override;

  /// Gets the cutoff wavenumber for a particular waveguide mode
  /// \param modeType The mode type to get (TE, TM, TEM)
  /// \param n The angular number of the mode
  /// \param m The radial number of the mode
  /// \Returns The cutoff wavenumber in m^-1
  double GetCutoffWavenumber(Mode_t modeType, unsigned int n,
                             unsigned int m) override;

  /// Gets the resonant frequency for a particle mode
  /// \param modeType The type of mode (TE or TM)
  /// \param n The angular number of the mode
  /// \param m The radial number of the mode
  /// \param l The mode number in the z direction of the waveguide
  /// \Returns The resonant frequency of the chosen mode (in Hertz)
  double GetResonantModeFrequency(Mode_t modeType, int n, int m,
                                  int l) override;

  /// Calculates the required normalisation of the electric fields
  /// If you multiply the fields by the result should give correct normalisation
  /// \param modeType The type of mode (TE or TM)
  /// \param n The angular number of the mode
  /// \param m The radial number of the mode
  /// \param omega Angular frequency for the chosen wave
  /// \param A The constant factor for one polarisation state
  /// \param B The constant factor for the other polarisation state
  /// \param nSurfPnts Number of points in each dimension to test
  double GetEFieldNormalisation(Mode_t modeType, int n, int m, double omega,
                                double A, double B, int nSurfPnts);

  double GetHFieldIntegral(Mode_t modeType, int n, int m, double omega,
                           double A, double B, int nSurfPnts);

  /// @brief Gets the field amplitude from a moving electron in the guide
  /// @param modeType The mode type to get (TE, TM)
  /// @param n The angular number of the mode
  /// @param m The radial number of the mode
  /// @param omega Angular frequency of the chosen wave
  /// @param ePos The electron position vector
  /// @param eVel The electron velocity vector
  /// @param normA Normalisation of one polarisation
  /// @param normB Normalisation of one polarisation
  /// @param isPositive Do we want the positive or negative amplitude?
  /// @return The field amplitude at a given time
  std::complex<double> GetFieldAmp(Mode_t modeType, unsigned int n,
                                   unsigned int m, double omega, TVector3 ePos,
                                   TVector3 eVel, double normA, double normB,
                                   bool isPositive) override;

  /// @brief Calculate and set Pn for a given circular waveguide mode
  /// @param modeType TE or TM mode
  /// @param n Angular mode number
  /// @param m Radial mode number
  /// @param omega Angular frequency of chosen mode
  /// @param A Constant 1
  /// @param B Constant 2
  /// @param nSurfPnts Number of points to use for integration
  void CalculatePn(Mode_t modeType, unsigned int n, unsigned int m,
                   double omega, unsigned int nSurfPnts, double A, double B);
};
}  // namespace rad

#endif
