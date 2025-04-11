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
  CircularWaveguide(double radius, double length);

  /// Return the inner radius of the waveguide
  /// \Returns The waveguide inner radius (in metres)
  double GetInnerRadius() { return a; }

  /// Return the length of the waveguide
  /// \Returns The waveguide length (in metres)
  double GetLength() { return d; }

  /// @brief Gets the electric field vector for a given mode at a point
  /// @param pos Position vector (in metres)
  /// @param mode The mode to get
  /// @param A Normalisation constant
  /// @param omega Angular frequency of the wave
  /// @param state Choose polarisation
  /// @return Electric field vector at that point
  TVector3 GetModeEField(TVector3 pos, WaveguideMode mode, double A,
                         double omega, bool state) override;

  /// @brief Gets the H field vector for a given mode at a point
  /// @param pos The position vector (in metres)
  /// @param mode The mode to get
  /// @param A Arbitrary amplitude for part of solution
  /// @param omega Angular frequency of the chosen wave
  /// @param state Choose polarisation state
  /// @return The mode H field vector at the supplied point
  TVector3 GetModeHField(TVector3 pos, WaveguideMode mode, double A,
                         double omega, bool state) override;

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
  /// @param nSurfPnts Number of points to use for integration
  void CalculatePn(WaveguideMode mode, double omega, unsigned int nSurfPnts);

  bool MultiplePolarisations() { return true; }
};
}  // namespace rad

#endif
