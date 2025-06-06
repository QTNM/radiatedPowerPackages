/*
  ICavity.h

  Base class for resonant cavities

  S. Jones 30-06-2023
*/

#ifndef ICAVITY_H
#define ICAVITY_H

#include <boost/math/special_functions/bessel.hpp>

#include "BasicFunctions/ComplexVector3.h"
#include "TVector3.h"

namespace rad {
class ICavity {
 public:
  virtual ~ICavity() {}

  // Different types of modes
  enum Mode_t { kTE, kTM, kTEM };

  /// @brief Getter function for probe position
  /// @return 3-vector of probe position
  TVector3 GetProbePosition() { return probePos; }

  /// @brief Calculates the frequency of a resonant mode
  /// @param modeType TE, TM or TEM
  /// @param n First mode index
  /// @param m Second mode index
  /// @param l Axial mode index
  /// @return Frequency in Hertz
  virtual double GetResonantModeF(Mode_t modeType, unsigned int n,
                                  unsigned int m, unsigned int l) = 0;

  /// @brief Gets the cutoff frequency for a particular waveguide mode
  /// @param modeType The mode type to get (TE, TM, TEM)
  /// @param n The first mode index to get
  /// @param m The second mode index to get
  /// @return The cutoff frequency in Hertz
  virtual double GetCutoffFrequency(Mode_t modeType, int n, int m) = 0;

  /// @brief Get mode electric field with no time variation
  /// @param pos Position to evaluate field at
  /// @param modeType Select TE, TM, or TEM
  /// @param A Normalisation factor
  /// @param i First mode index
  /// @param j Second mode index
  /// @param k Third mode index
  /// @param state Choose from two polarisation states
  /// @return Complex vector of field
  virtual ComplexVector3 GetModalEField(TVector3 pos, Mode_t modeType, double A,
                                        unsigned int i, unsigned int j,
                                        unsigned int k, bool state) = 0;

  /// @brief Calculate factor to correctly normalise mode fields
  /// @param modeType Select TE or TM
  /// @param i First mode index
  /// @param j Second mode index
  /// @param k Third mode index
  /// @param state Choose polarisation state
  /// @return Factor which unnormalised field must be multiplied by to give
  /// correct integral
  virtual std::complex<double> GetModeNormalisation(Mode_t modeType,
                                                    unsigned int i,
                                                    unsigned int j,
                                                    unsigned int k,
                                                    bool state) = 0;

 private:
  TVector3 probePos;

 protected:
  /// @brief Probe position setter
  /// @param probe 3-vector of desired probe position
  void SetProbePosition(TVector3 probe) { probePos = probe; }
};
}  // namespace rad

#endif