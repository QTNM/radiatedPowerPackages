/*
  Circular resonant cavity

  S. Jones 30-06-2023
*/

#ifndef CIRCULAR_CAVITY_H
#define CIRCULAR_CAVITY_H

#include "Waveguides/ICavity.h"

namespace rad {
class CircularCavity : public ICavity {
 private:
  double a;
  double l;

 public:
  /// Parametrised constructor
  /// @param radius Cavity radius in metres
  /// @param length Cavity length in metres
  CircularCavity(double radius, double length) : a(radius), l(length) {}

  double GetResonantModeF(Mode_t modeType, unsigned int n, unsigned int m,
                          unsigned int l) override;
};
}  // namespace rad

#endif