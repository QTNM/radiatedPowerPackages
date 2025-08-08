// BasicFunctions.h - Consolidated header for all utility functions
#ifndef BASIC_FUNCTIONS_H
#define BASIC_FUNCTIONS_H

// Include ROOT-based utilities (requires ROOT)
#include "utilities/ROOTUtils/CorrelationAnalysis.h"
#include "utilities/ROOTUtils/FFTAnalysis.h"
#include "utilities/ROOTUtils/GraphUtils.h"
#include "utilities/ROOTUtils/HistogramUtils.h"

// Legacy includes for backward compatibility
#include "Math/Point3D.h"
#include "Math/Vector3D.h"
#include "TGraph.h"
#include "TH1.h"
#include "TH2.h"
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
