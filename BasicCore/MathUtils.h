// MathUtils.h
// Simple mathematical utility functions - header-only with no external
// dependencies

#ifndef MATHUTILS_H
#define MATHUTILS_H

#include "Constants.h"

namespace rad {

/// @brief Convert degrees to radians
/// @param degrees Angle in degrees
/// @return Angle in radians
inline constexpr double DegToRad(double degrees) {
  return degrees * PI / 180.0;
}

/// @brief Convert radians to degrees
/// @param radians Angle in radians
/// @return Angle in degrees
inline constexpr double RadToDeg(double radians) {
  return radians * 180.0 / PI;
}

}  // namespace rad

#endif