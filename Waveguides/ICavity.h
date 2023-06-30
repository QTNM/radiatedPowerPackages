/*
  ICavity.h

  Base class for resonant cavities

  S. Jones 30-06-2023
*/

#ifndef ICAVITY_H
#define ICAVITY_H

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
};
}  // namespace rad

#endif