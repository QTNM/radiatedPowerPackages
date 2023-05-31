/*
  BaseField.h

  Abstract base class for magnetic fields used electron tracking.
  Based on transcription of T. Goffrey's original python code
*/

#ifndef BASE_FIELD_H
#define BASE_FIELD_H

#include "TVector3.h"

namespace rad {
class BaseField {
 public:
  /// @brief Virtual function for calculating electric field
  /// @param v Point at which to evaluate electric field
  /// @return Electric field vector in volts/metre
  virtual TVector3 evaluate_e_field_at_point(TVector3 v) = 0;

  /// @brief Calculates electric field magnitude at a point
  /// @param v Point at which to evaluate field
  /// @return Electric field magnitude in volts/metre
  double evaluate_e_field_magnitude(TVector3 v);

  /// @brief Virtual function for calculating magnetic field
  /// @param v Point at which to evaluate magnetic field
  /// @return Magnetic field vector in tesla
  virtual TVector3 evaluate_field_at_point(const TVector3 v) = 0;

  /// @brief Magnetic field vector magnitude
  /// @param v Point at which to evaluate field
  /// @return Magnetic field magnitude in Tesla
  double evaluate_field_magnitude(TVector3 v);

  virtual ~BaseField() {}
};

}  // namespace rad

#endif
