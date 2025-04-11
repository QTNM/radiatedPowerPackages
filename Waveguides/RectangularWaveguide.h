/*
  RectangularWaveguide.h

  Derived class for a rectangular waveguide
*/

#ifndef RECTANGULAR_WAVEGUIDE_H
#define RECTANGULAR_WAVEGUIDE_H

#include "Waveguides/IWaveguide.h"

namespace rad {

/// Class describing a rectangular waveguide
/// Centre of waveguide located at x = y = z = 0
class RectangularWaveguide : public IWaveguide {
 private:
  // Convention is that a > b and that longest side of waveguide is x-axis
  double a;  // Side length in metres
  double b;  // Side length in metres
  double d;  // Axial length in metres

 public:
  /// @brief Parametrised constructor
  /// @param longSide The longer side of the waveguide (in metres)
  /// @param shortSide The shorter side of the waveguide (in metres)
  /// @param length The axial length of the waveguide (in metres);
  RectangularWaveguide(double longSide, double shortSide, double length);

  /// Gets the long dimension of the waveguide
  /// \Returns The long (x) dimension of the waveguide (in metres)
  double GetLongDimension() { return a; }

  /// Gets the short dimension of the waveguide
  /// \Returns The short (y) dimension of the waveguide (in metres)
  double GetShortDimension() { return b; }

  /// @brief Gets the real electric field vector for a given mode at a point
  /// @param pos The position vector (in metres)
  /// @param mode The mode type to get (either TE or TM)
  /// @param A Arbitrary amplitude for solution
  /// @param omega Angular frequency of the chosen wave
  /// @param state Does nothing in this case
  /// @return The mode electric field vector at the supplied point
  TVector3 GetModeEField(TVector3 pos, WaveguideMode mode, double A,
                         double omega, bool state) override;

  /// @brief Gets the H field vector for a given mode at a point
  /// @param pos The position vector (in metres)
  /// @param modeType The mode type to get (either TE or TM)
  /// @param A Arbitrary amplitude for solution
  /// @param omega Angular frequency of the chosen wave
  /// @param state Choose polarisation state (not applicable here)
  /// @return The mode H field vector at the supplied point
  TVector3 GetModeHField(TVector3 pos, WaveguideMode mode, double A,
                         double omega, bool state) override;

  /// Gets the cutoff frequency for a particular mode
  /// \param modeType The mode type to use
  /// \Returns The cutoff frequency of the mode in Hertz
  double GetCutoffFrequency(WaveguideMode mode) override;

  /// @brief Gets the cutoff wavenumber for a given mode
  /// @param mode The mode type to use
  /// @return The cutoff wavenumber in units of m^{-1}
  double GetCutoffWavenumber(WaveguideMode mode) override;

  /// @brief Calculates the normalisation integral of the mode electric field
  /// @param mode The mode (TE or TM)
  /// @param omega Angular frequency for the chosen wave
  /// @param A The constant factor
  /// @param nSurfPnts Number of points in each dimension to test
  /// @return Electric field integral
  double GetEFieldIntegral(WaveguideMode mode, double omega, double A,
                           int nSurfPnts, bool state) override;

  /// @brief Gets the field amplitude from a moving electron in the guide
  /// @param mode The mode to do the calculation for
  /// @param omega Angular frequency of the chosen wave
  /// @param ePos The electron position vector
  /// @param eVel The electron velocity vector
  /// @param normA Normalisation
  /// @param state Does nothing here
  /// @param isPositive Do we want the positive or negative amplitude?
  /// @return The field amplitude at a given time
  double GetFieldAmp(WaveguideMode mode, double omega, TVector3 ePos,
                     TVector3 eVel, double normA, bool state,
                     bool isPositive) override;

  /// @brief Calculate and set Pn for a given mode
  /// @param modeType TE or TM mode
  /// @param omega Angular frequency
  /// @param nSurfPnts Number of points for integration
  void CalculatePn(WaveguideMode mode, double omega, unsigned int nSurfPnts);

  bool MultiplePolarisations() { return false; }
};
}  // namespace rad

#endif
