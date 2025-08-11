// BasicFunctions.h
#ifndef BASIC_FUNCTIONS_H
#define BASIC_FUNCTIONS_H

#include "TVector3.h"

namespace rad {

/// @brief Rotate vector to global coordinate system
/// @param v The vector to be rotated
/// @param xAx The X axis of the frame you are rotating from
/// @param yAx The Y axis of the frame you are rotating from
/// @param zAx The Z axis of the frame you are rotating from
/// @return The rotated vector
TVector3 RotateToGlobalCoords(TVector3 v, TVector3 xAx, TVector3 yAx,
                              TVector3 zAx);

/// @brief Rotate vector to different coordinate system
/// @param v Vector to be rotated
/// @param newX X axis of frame to be rotated to
/// @param newY Y axis of frame to be rotated to
/// @param newZ Z axis of frame to be rotated to
/// @return 3-vector of rotated vector
TVector3 RotateToCoords(TVector3 v, TVector3 newX, TVector3 newY,
                        TVector3 newZ);
}  // namespace rad

#endif
