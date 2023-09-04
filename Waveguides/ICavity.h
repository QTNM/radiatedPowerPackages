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

 private:
  TVector3 probePos;

 protected:
  /// @brief Getter function for probe position
  /// @return 3-vector of probe position
  TVector3 GetProbePosition() { return probePos; }

  /// @brief Probe position setter
  /// @param probe 3-vector of desired probe position
  void SetProbePosition(TVector3 probe) { probePos = probe; }
};
}  // namespace rad

#endif