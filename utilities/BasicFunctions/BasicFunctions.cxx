/// BasicFunctions.cxx - Non-ROOT utility functions
#include "utilities/BasicFunctions/BasicFunctions.h"

#include "utilities/BasicCore/Constants.h"
#include "TVector3.h"

TVector3 rad::RotateToGlobalCoords(TVector3 v, TVector3 xAx, TVector3 yAx,
                                   TVector3 zAx) {
  // We are just transforming to the ROOT frame so our new axis coordinates
  // are just the unit vectors
  TVector3 unitX(1, 0, 0);
  TVector3 unitY(0, 1, 0);
  TVector3 unitZ(0, 0, 1);
  double x1Prime{unitX.Dot(xAx.Unit()) * v.X() + unitX.Dot(yAx.Unit()) * v.Y() +
                 unitX.Dot(zAx.Unit()) * v.Z()};
  double x2Prime{unitY.Dot(xAx.Unit()) * v.X() + unitY.Dot(yAx.Unit()) * v.Y() +
                 unitY.Dot(zAx.Unit()) * v.Z()};
  double x3Prime{unitZ.Dot(xAx.Unit()) * v.X() + unitZ.Dot(yAx.Unit()) * v.Y() +
                 unitZ.Dot(zAx.Unit()) * v.Z()};
  TVector3 newVector{x1Prime * unitX + x2Prime * unitY + x3Prime * unitZ};
  return newVector;
}

TVector3 rad::RotateToCoords(TVector3 v, TVector3 newX, TVector3 newY,
                             TVector3 newZ) {
  // We are just transforming from the ROOT frame so our old axis coordinates
  // are just the unit vectors
  TVector3 oldX(1, 0, 0);
  TVector3 oldY(0, 1, 0);
  TVector3 oldZ(0, 0, 1);
  double pXPrime{newX.X() * v.X() + newY.X() * v.Y() + newZ.X() * v.Z()};
  double pYPrime{newX.Y() * v.X() + newY.Y() * v.Y() + newZ.Y() * v.Z()};
  double pZPrime{newX.Z() * v.X() + newY.Z() * v.Y() + newZ.Z() * v.Z()};
  TVector3 newVector(pXPrime, pYPrime, pZPrime);
  newVector = newVector.Unit();
  return newVector;
}