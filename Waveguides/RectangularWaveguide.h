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
  /// Parametrised constructor
  /// \param longSide The longer side of the waveguide (in metres)
  /// \param shortSide The shorter side of the waveguide (in metres)
  /// \param length The axial length of the waveguide (in metres);
  RectangularWaveguide(double longSide, double shortSide, double length);

  /// Gets the long dimension of the waveguide
  /// \Returns The long (x) dimension of the waveguide (in metres)
  double GetLongDimension() { return a; }

  /// Gets the short dimension of the waveguide
  /// \Returns The short (y) dimension of the waveguide (in metres)
  double GetShortDimension() { return b; }

  /// Gets the complex electric field vector for a given mode at a point
  /// \param modeType The mode type to get (either TE or TM)
  /// \param m The mode number in the x direction of the waveguide
  /// \param n The mode number in the y direction of the waveguide
  /// \param pos The position vector (in metres)
  /// \param omega Angular frequency of the chosen wave
  /// \param A Arbitrary amplitude for solution (default = 1)
  /// \param B Arbitrary amplitude for solution has no effect
  /// \Returns The mode electric field vector at the supplied point
  ComplexVector3 GetModeEFieldComplex(Mode_t modeType, int m, int n,
                                      TVector3 pos, double omega, double A = 1,
                                      double B = 0) override;

  /// Gets the real electric field vector for a given mode at a point
  /// \param modeType The mode type to get (either TE or TM)
  /// \param m The mode number in the x direction of the waveguide
  /// \param n The mode number in the y direction of the waveguide
  /// \param pos The position vector (in metres)
  /// \param omega Angular frequency of the chosen wave
  /// \param A Arbitrary amplitude for solution (default = 1)
  /// \param B Arbitrary amplitude for solution has no effect
  /// \Returns The mode electric field vector at the supplied point
  TVector3 GetModeEField(Mode_t modeType, int m, int n, TVector3 pos,
                         double omega, double A = 1, double B = 0) override;

  TVector3 GetModalEField(Mode_t modeType, int m, int n, TVector3 pos,
                          double omega, double A = 1);

  ComplexVector3 GetNormalisedEField(Mode_t modeType, int m, int n,
                                     TVector3 pos, double omega);

  /// Gets the complex magnetic field strength vector for a given mode at a
  /// point \param modeType The mode type to get (either TE or TM) \param m The
  /// mode number in the x direction of the waveguide \param n The mode number
  /// in the y direction of the waveguide \param pos The position vector (in
  /// metres) \param omega Angular frequency of the chosen wave \param A
  /// Arbitrary amplitude for solution (default = 1) \Returns The mode H field
  /// vector at the supplied point
  ComplexVector3 GetModeHFieldComplex(Mode_t modeType, int m, int n,
                                      TVector3 pos, double omega, double A = 1,
                                      double B = 0) override;

  /// @brief Gets the H field vector for a given mode at a point
  /// @param modeType The mode type to get (either TE or TM)
  /// @param m The mode number in the x direction of the waveguide
  /// @param n The mode number in the y direction of the waveguide
  /// @param pos The position vector (in metres)
  /// @param omega Angular frequency of the chosen wave
  /// @param A Arbitrary amplitude for solution (default = 1)
  /// @return The mode H field vector at the supplied point
  TVector3 GetModeHField(Mode_t modeType, int m, int n, TVector3 pos,
                         double omega, double A = 1, double B = 0) override;

  TVector3 GetModalHField(Mode_t modeType, int m, int n, TVector3 pos,
                          double omega, double A = 1);

  ComplexVector3 GetNormalisedHField(Mode_t modeType, int m, int n,
                                     TVector3 pos, double omega);

  /// Gets the cutoff frequency for a particular mode
  /// \param modeType The mode type to use (either TE or TM)
  /// \param m The mode number in the x direction of the waveguide
  /// \param n The mode number in the y direction of the waveguide
  /// \Returns The cutoff frequency of the mode in Hertz
  double GetCutoffFrequency(Mode_t modeType, int m, int n) override;

  /// @brief Gets the cutoff wavenumber for a given mode
  /// @param modeType The mode type to use (this doesn't matter for rectangular
  /// guides)
  /// @param m Mode order in the x direction of the waveguide
  /// @param n Mode order in the y direction of the waveguide
  /// @return The cutoff wavenumber in units of m^{-1}
  double GetCutoffWavenumber(Mode_t modeType, unsigned int m,
                             unsigned int n) override;

  /// Gets the resonant frequency for a particle mode
  /// \param modeType The type of mode (TE or TM)
  /// \param m The mode number in the x direction of the waveguide
  /// \param n The mode number in the y direction of the waveguide
  /// \param l The mode number in the z direction of the waveguide
  /// \Returns The resonant frequency of the chosen mode (in Hertz)
  double GetResonantModeFrequency(Mode_t modeType, int m, int n,
                                  int l) override;

  /// Gets the field amplitude from a moving electron in the guide
  /// \param modeType The mode type to get (TE, TM)
  /// \param m The mode number in the x direction of the waveguide
  /// \param n The mode number in the y direction of the waveguide
  /// \param omega Angular frequency of the chosen wave
  /// \param ePos The electron position vector
  /// \param eVel The electron velocity vector
  /// \param normA Normalisation of one polarisation (not used here)
  /// \param normB Normalisation of one polarisation (not used here)
  /// \Returns The field amplitude at a given time
  std::complex<double> GetPositiveFieldAmp(Mode_t modeType, unsigned int m,
                                           unsigned int n, double omega,
                                           TVector3 ePos, TVector3 eVel,
                                           double normA, double normB) override;

  /// @brief Calculate and set Pn for a given mode
  /// @param modeType TE or TM mode
  /// @param n
  /// @param m
  /// @param omega Angular frequency
  /// @param nSurfPnts Number of points for integration
  void CalculatePn(Mode_t modeType, unsigned int n, unsigned int m,
                   double omega, unsigned int nSurfPnts);
};
}  // namespace rad

#endif
