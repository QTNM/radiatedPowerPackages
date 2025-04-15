/*
  IWaveguide.h

  Abstract base class for waveguides
*/

#ifndef IWAVEGUIDE_H
#define IWAVEGUIDE_H

#include "BasicFunctions/ComplexVector3.h"
#include "TMath.h"
#include "TVector3.h"
#include "Waveguides/WaveguideMode.h"

namespace rad {

class IWaveguide {
 public:
  IWaveguide() : Pn(0) {}

  virtual ~IWaveguide() {}

  // Different kinds of modes we may have
  enum Mode_t { kTE, kTM, kTEM };

  /// @brief Gets the electric field vector for a given mode at a point
  /// @param pos The position vector (in metres)
  /// @param mode The waveguide mode to get
  /// @param A Arbitrary amplitude for solution
  /// @param omega Angular frequency of the chosen wave
  /// @param state Choose polarisation state (where applicable)
  /// @return The mode electric field vector at the supplied point
  virtual TVector3 GetModeEField(TVector3 pos, WaveguideMode mode, double A,
                                 double omega, bool state) = 0;

  /// @brief Gets the H field vector for a given mode at a point
  /// @param pos The position vector (in metres)
  /// @param mode The mode to get
  /// @param A Arbitrary amplitude for part of solution
  /// @param omega Angular frequency of the chosen wave
  /// @param state Choose polarisation state (where applicable)
  /// @return The mode H field vector at the supplied point
  virtual TVector3 GetModeHField(TVector3 pos, WaveguideMode mode, double A,
                                 double omega, bool state) = 0;

  /// @brief Gets the characteristic mode impedance for a given mode
  /// @param mode The mode to get
  /// @param omega Angular frequency of the chosen wave
  /// @return The impedance of the mode (in Ohms)
  double GetModeImpedance(WaveguideMode mode, double omega);

  /// @brief Calculates the cutoff frequency for a particular waveguide mode
  /// @param mode The mode to do the calculation for
  /// @return The cutoff frequency in Hertz
  virtual double GetCutoffFrequency(WaveguideMode mode) = 0;

  /// @brief Calculates the cutoff wavenumber for a particular waveguide mode
  /// @param mode The mode to do the calculation for
  /// @return The cutoff wavenumber in m^-1
  virtual double GetCutoffWavenumber(WaveguideMode mode) = 0;

  /// @brief Check if a mode propagates at a given frequency
  /// @param mode The mode to be checked
  /// @param f Frequency of wave in guide in Hertz
  /// @return True is mode propagates, false if not
  bool ModePropagates(WaveguideMode mode, double f);

  /// @brief Gets the field amplitude from a moving electron in the guide
  /// @param modeType The mode to calculate for
  /// @param omega Angular frequency of the chosen wave
  /// @param ePos The electron position vector
  /// @param eVel The electron velocity vector
  /// @param norm Constant factor multiplying mode
  /// @param state Choose polarisation state (where more than one exists)
  /// @param isPositive Do we want the positive or negative amplitude?
  /// @return The field amplitude at a given time
  virtual double GetFieldAmp(WaveguideMode mode, double omega, TVector3 ePos,
                             TVector3 eVel, double norm, bool state,
                             bool isPositive) = 0;

  /// @brief Getter function for normalisation constant
  /// @return Normalisation constant
  double GetPn() { return Pn; }

  /// @brief Calculates the required normalisation of the electric fields. If
  /// you multiply the fields by the result should give correct normalisation
  /// @param mode The mode
  /// @param omega Angular frequency for the chosen wave
  /// @param A The constant factor
  /// @param nSurfPnts Number of points in each dimension to test
  /// @param state Can be used for choosing polarisation state
  virtual double GetEFieldIntegral(WaveguideMode mode, double omega, double A,
                                   int nSurfPnts, bool state) = 0;

  /// @brief Function for determining if we need to calculate multiple
  /// polarisations
  /// @return True is there are multiple polarisations
  virtual bool MultiplePolarisations() = 0;

  /// @brief Function for calculating the phase velocity of a mode
  /// @param mode The waveguide mode to calculate for
  /// @param f The frequency for which to calculate the phase velocity in Hertz
  /// @return The phase velocity in m/s
  double GetPhaseVelocity(WaveguideMode mode, double f);

  /// @brief Function for calculating the group velocity of a mode
  /// @param mode The waveguide mode to calculate for
  /// @param f The frequency for which to calculate the group velocity in Hertz
  /// @return The group velocity in m/s
  double GetGroupVelocity(WaveguideMode mode, double f);

  /// @brief Function for calculating normalisation constant Pn for a given
  /// mode. Also sets the variable Pn for the waveguide
  /// @param mode Mode type in questions
  /// @param omega Angular frequency of the chosen wave
  /// @param nSurfPnts Number of points in each dimension to test
  virtual void CalculatePn(WaveguideMode mode, double omega,
                           unsigned int nSurfPnts = 100) = 0;

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
