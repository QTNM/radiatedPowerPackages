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

 public:
  /// @brief Parametrised  constructor
  /// @param radius Inner radius of the waveguide (in metres)
  /// @param length Length of the waveguide (in metres)
  /// @param probePos Probe position 3-vector
  CircularWaveguide(double radius, double length, TVector3 probePos);

  /// Return the inner radius of the waveguide
  /// \Returns The waveguide inner radius (in metres)
  double GetInnerRadius() { return a; }

  /// Return the length of the waveguide
  /// \Returns The waveguide length (in metres)
  double GetLength() { return d; }

  /// Complex electric field vector
  /// \param mode The mode to get
  /// \param pos The position vector (in metres)
  /// \param omega Angular frequency of the chosen wave
  /// \param A Arbitrary amplitude for part of solution (default = 1)
  /// \param B Arbitrary amplitude for part of solution (default = 0)
  /// \Returns The mode complex electric field vector at the supplied point
  ComplexVector3 GetModeEFieldComplex(WaveguideMode mode, TVector3 pos,
                                      double omega, double A = 1,
                                      double B = 0) override;

  /// @brief Gets the electric field vector for a given mode at a point
  /// @param pos Position vector (in metres)
  /// @param mode The mode to get
  /// @param A Normalisation constant
  /// @param omega Angular frequency of the wave
  /// @param state Choose polarisation
  /// @return Electric field vector at that point
  TVector3 GetModeEField(TVector3 pos, WaveguideMode mode, double A,
                         double omega, bool state) override;

  /// @brief Gets the complex magnetic field strength vector for a given mode at
  /// a point
  /// @param mode The mode to get
  /// @param pos The position vector (in metres)
  /// @param omega Angular frequency of the chosen wave
  /// @param A Arbitrary amplitude for part of solution (default = 1)
  /// @param B Arbitrary amplitude for part of solution (default = 0)
  /// @return The mode H field vector at the supplied point
  ComplexVector3 GetModeHFieldComplex(WaveguideMode mode, TVector3 pos,
                                      double omega, double A = 1,
                                      double B = 0) override;

  /// Gets the H field vector for a given mode at a point
  /// \param mode The mode to get
  /// \param pos The position vector (in metres)
  /// \param omega Angular frequency of the chosen wave
  /// \param A Arbitrary amplitude for part of solution (default = 1)
  /// \param B Arbitrary amplitude for part of solution (default = 0)
  /// \Returns The mode H field vector at the supplied point
  TVector3 GetModeHField(WaveguideMode mode, TVector3 pos, double omega,
                         double A = 1, double B = 0) override;

  ComplexVector3 GetModalHField(WaveguideMode mode, TVector3 pos, double omega,
                                double A = 1, double B = 0);

  /// @brief Gets the cutoff frequency for a particular waveguide mode
  /// \param mode The mode to get
  /// @return The cutoff frequency in Hertz
  double GetCutoffFrequency(WaveguideMode mode) override;

  /// @brief Gets the cutoff wavenumber for a particular waveguide mode
  /// @param modeType The mode to get
  /// @return The cutoff wavenumber in m^-1
  double GetCutoffWavenumber(WaveguideMode mode) override;

  /// @brief Calculates the normalisation integral of the mode electric field
  /// @param modeType The type of mode
  /// @param omega Angular frequency for the chosen wave
  /// @param A The constant factor
  /// @param nSurfPnts Number of points in each dimension to test
  /// @return Electric field integral
  double GetEFieldIntegral(WaveguideMode mode, double omega, double A,
                           int nSurfPnts, bool state) override;

  double GetHFieldIntegral(WaveguideMode mode, double omega, double A, double B,
                           int nSurfPnts);

  /// @brief Gets the field amplitude from a moving electron in the guide
  /// @param mode The mode to calculate for
  /// @param omega Angular frequency of the chosen wave
  /// @param ePos The electron position vector
  /// @param eVel The electron velocity vector
  /// @param normA Normalisation
  /// @param state Choose polarisation state
  /// @param isPositive Do we want the positive or negative amplitude?
  /// @return The field amplitude at a given time
  double GetFieldAmp(WaveguideMode mode, double omega, TVector3 ePos,
                     TVector3 eVel, double normA, bool state,
                     bool isPositive) override;

  /// @brief Calculate and set Pn for a given circular waveguide mode
  /// @param mode The mode
  /// @param omega Angular frequency of chosen mode
  /// @param A Constant 1
  /// @param B Constant 2
  /// @param nSurfPnts Number of points to use for integration
  void CalculatePn(WaveguideMode mode, double omega, unsigned int nSurfPnts,
                   double A, double B);

  bool MultiplePolarisations() { return true; }
};
}  // namespace rad

#endif
