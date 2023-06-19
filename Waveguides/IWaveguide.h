/*
  IWaveguide.h

  Abstract base class for waveguides
*/

#ifndef IWAVEGUIDE_H
#define IWAVEGUIDE_H

#include "BasicFunctions/ComplexVector3.h"
#include "TMath.h"
#include "TVector3.h"

namespace rad {

class IWaveguide {
 public:
  IWaveguide() : Pn(0) {}

  virtual ~IWaveguide() {}

  // Different kinds of modes we may have
  enum Mode_t { kTE, kTM, kTEM };

  /// Gets the complex electric field vector for a given mode at a point
  /// \param modeType The mode type to get (TE, TM, TEM)
  /// \param pos The position vector (in metres)
  /// \param omega Angular frequency of the chosen wave
  /// \param A Arbitrary amplitude for part of solution (default = 1)
  /// \param B Arbitrary amplitude for part of solution (default = 0)
  /// \Returns The mode electric field vector at the supplied point
  virtual ComplexVector3 GetModeEFieldComplex(Mode_t modeType, int n, int m,
                                              TVector3 pos, double omega,
                                              double A = 1, double B = 0) = 0;

  /// Gets the electric field vector for a given mode at a point
  /// \param modeType The mode type to get (TE, TM, TEM)
  /// \param n The angular number of the mode
  /// \param m The radial number of the mode
  /// \param pos The position vector (in metres)
  /// \param omega Angular frequency of the chosen wave
  /// \param A Arbitrary amplitude for part of solution (default = 1)
  /// \param B Arbitrary amplitude for part of solution (default = 0)
  /// \Returns The mode electric field vector at the supplied point
  virtual TVector3 GetModeEField(Mode_t modeType, int n, int m, TVector3 pos,
                                 double omega, double A = 1, double B = 0) = 0;

  /// Gets the complex magnetic field strength vector for a given mode at a
  /// point \param modeType The mode type to get (TE, TM, TEM) \param n The
  /// angular number of the mode \param m The radial number of the mode \param
  /// pos The position vector (in metres) \param omega Angular frequency of the
  /// chosen wave \param A Arbitrary amplitude for part of solution (default =
  /// 1) \param B Arbitrary amplitude for part of solution (default = 0)
  /// \Returns The mode H field vector at the supplied point
  virtual ComplexVector3 GetModeHFieldComplex(Mode_t modeType, int n, int m,
                                              TVector3 pos, double omega,
                                              double A = 1, double B = 0) = 0;

  /// Gets the H field vector for a given mode at a point
  /// \param modeType The mode type to get (TE, TM, TEM)
  /// \param n The angular number of the mode
  /// \param m The radial number of the mode
  /// \param pos The position vector (in metres)
  /// \param omega Angular frequency of the chosen wave
  /// \param A Arbitrary amplitude for part of solution (default = 1)
  /// \param B Arbitrary amplitude for part of solution (default = 0)
  /// \Returns The mode H field vector at the supplied point
  virtual TVector3 GetModeHField(Mode_t modeType, int n, int m, TVector3 pos,
                                 double omega, double A = 1, double B = 0) = 0;

  /// Gets the characteristic mode impedance for a given mode
  /// \param modeType The mode type to get (TE, TM, TEM)
  /// \param n The first mode index
  /// \param m The second mode index
  /// \param omega Angular frequency of the chosen wave
  /// \Returns The impedance of the mode (in Ohms)
  double GetModeImpedance(Mode_t modeType, unsigned int n, unsigned int m,
                          double omega);

  /// Gets the cutoff frequency for a particular waveguide mode
  /// \param modeType The mode type to get (TE, TM, TEM)
  /// \param n The first mode index to get
  /// \param m The second mode index to get
  /// \Returns The cutoff frequency in Hertz
  virtual double GetCutoffFrequency(Mode_t modeType, int n, int m) = 0;

  /// Gets the cutoff wavenumber for a particular waveguide mode
  /// \param modeType The mode type to get (TE, TM, TEM)
  /// \param n The first mode index to get
  /// \param m The second mode index to get
  /// \Returns The cutoff wavenumber in m^-1
  virtual double GetCutoffWavenumber(Mode_t modeType, unsigned int n,
                                     unsigned int m) = 0;

  /// Gets the frequency of a particular resonant mode
  /// \Returns The resonant frequency of the chosen mode (in Hertz)
  virtual double GetResonantModeFrequency(Mode_t modeType, int n, int m,
                                          int l) = 0;

  /// Gets the field amplitude from a moving electron in the guide
  /// \param modeType The mode type to get (TE, TM, TEM)
  /// \param n The first mode index to get
  /// \param m The second mode index to get
  /// \param omega Angular frequency of the chosen wave
  /// \param ePos The electron position vector
  /// \param eVel The electron velocity vector
  /// \param normA Normalisation of one polarisation (circular guide only)
  /// \param normB Normalisation of one polarisation (circular guide only)
  /// \Returns The field amplitude at a given time
  virtual std::complex<double> GetPositiveFieldAmp(
      Mode_t modeType, unsigned int n, unsigned int m, double omega,
      TVector3 ePos, TVector3 eVel, double normA, double normB) = 0;

  /// @brief Calculate a normalisation constant
  /// @param modeType TE or TM
  /// @param n
  /// @param m
  /// @param omega Angular frequency of radiation
  /// @param nSurfPnts Number of points to use for integration
  /// @return Value of Pn
  virtual void CalculatePn(Mode_t modeType, unsigned int n, unsigned int m,
                           double omega, unsigned int nSurfPnts) = 0;

  /// @brief Getter function for normalisation constant
  /// @return Normalisation constant
  double GetPn() { return Pn; }

 private:
  // Normalisation constant
  double Pn;

 protected:
  /// @brief Setter function for normalisation constant
  /// @param val
  void SetPn(double val) { Pn = val; }
};
}  // namespace rad

#endif
