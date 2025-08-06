// Physics.h
// Header-only inline functions for basic physics calculations

#ifndef PHYSICS_H
#define PHYSICS_H

#include <cmath>

#include "Constants.h"

namespace rad {
/// @brief Calculate retarded time for electromagnetic field calculations
/// @param fieldPointPos Field point position (x, y, z) in meters
/// @param ePositionPos Source position (x, y, z) in meters
/// @param labTime Lab time in seconds
/// @return Retarded time in seconds
inline double CalcRetardedTime(const double fieldPointPos[3],
                               const double ePositionPos[3],
                               const double labTime) {
  double dx = ePositionPos[0] - fieldPointPos[0];
  double dy = ePositionPos[1] - fieldPointPos[1];
  double dz = ePositionPos[2] - fieldPointPos[2];
  double distance = std::sqrt(dx * dx + dy * dy + dz * dz);
  return labTime - distance / C;
}

/// @brief Calculate lab time from retarded time
/// @param fieldPointPos Field point position (x, y, z) in meters
/// @param ePositionPos Source position (x, y, z) in meters
/// @param tRet Retarded time in seconds
/// @return Lab time in seconds
inline double CalcTimeFromRetardedTime(const double fieldPointPos[3],
                                       const double ePositionPos[3],
                                       double tRet) {
  double dx = ePositionPos[0] - fieldPointPos[0];
  double dy = ePositionPos[1] - fieldPointPos[1];
  double dz = ePositionPos[2] - fieldPointPos[2];
  double distance = std::sqrt(dx * dx + dy * dy + dz * dz);
  return tRet + distance / C;
}
}  // namespace rad

#endif