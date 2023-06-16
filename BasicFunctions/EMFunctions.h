/// EMFunctions.h
/// Functions used for calculations of electromagnetic fields
#ifndef EM_FUNCTIONS_H
#define EM_FUNCTIONS_H

#include "BasicFunctions/Constants.h"

#include "TVector3.h"
#include "Math/Vector3D.h"
#include "Math/Point3D.h"

namespace rad
{
  /// Calculates electric field from a moving electron at a point
  /// Coordinates are all in same reference framce
  /// \param fieldPoint Vector of field point
  /// \param ePosition Vector of electron position
  /// \param eVelocity Electron velocity vector
  /// \param eAcceleration Electron acceleration vector
  /// \return The electric field vector
  ROOT::Math::XYZVector CalcEField(const ROOT::Math::XYZPoint fieldPoint, const ROOT::Math::XYZPoint ePosition, const ROOT::Math::XYZVector eVelocity, const ROOT::Math::XYZVector eAcceleration);
  /// Calculates electric field using TVector3 framework
  ROOT::Math::XYZVector CalcEField(const TVector3 fieldPoint, const TVector3 ePosition, const TVector3 eVelocity, const TVector3 eAcceleration);
  /// Calculate the far field component only
  TVector3 CalcEFarField(const TVector3 fieldPoint, const TVector3 ePosition, const TVector3 eVelocity, const TVector3 eAcceleration);
  
  /// \brief Calculates the non-relativistic electric field vector
  /// \param ePos Vector of electron position
  /// \param eVel Electron velocity vector
  /// \param eAcc Electron acceleration vector
  /// \return The electric field vector in V/m
  TVector3 CalcEFieldNR(TVector3 fp, TVector3 ePos, TVector3 eVel, TVector3 eAcc);

  /// Calculates magnetic field from a moving electron at a point
  /// Coordinates are all in same reference framce
  /// \param fieldPoint Vector of field point
  /// \param ePosition Vector of electron position
  /// \param eVelocity Electron velocity vector
  /// \param eAcceleration Electron acceleration vector
  /// \Returns The magnetic field vector
  ROOT::Math::XYZVector CalcBField(const ROOT::Math::XYZPoint fieldPoint, const ROOT::Math::XYZPoint ePosition, const ROOT::Math::XYZVector eVelocity, const ROOT::Math::XYZVector eAcceleration);
  /// Calculates magnetic field using TVector3 framework
  ROOT::Math::XYZVector CalcBField(const TVector3 fieldPoint, const TVector3 ePosition, const TVector3 eVelocity, const TVector3 eAcceleration); 
  /// Calculate the far field component only
  TVector3 CalcBFarField(const TVector3 fieldPoint, const TVector3 ePosition, const TVector3 eVelocity, const TVector3 eAcceleration);
  
  /// \brief Calculates the non-relativistic magnetic field vector
  /// \param ePos Vector of electron position
  /// \param eVel Electron velocity vector
  /// \param eAcc Electron acceleration vector
  /// \return The magnetic field vector 
  TVector3 CalcBFieldNR(TVector3 fp, TVector3 ePos, TVector3 eVel, TVector3 eAcc);

  /// Calculates Poynting vector from a moving electron at a point
  /// Coordinates are all in same reference framce
  /// \param fieldPoint Vector of field point
  /// \param ePosition Vector of electron position
  /// \param eVelocity Electron velocity vector
  /// \param eAcceleration Electron acceleration vector
  /// \Returns The Poynting vector
  ROOT::Math::XYZVector CalcPoyntingVec(const ROOT::Math::XYZPoint fieldPoint, const ROOT::Math::XYZPoint ePosition, const ROOT::Math::XYZVector eVelocity, const ROOT::Math::XYZVector eAcceleration);
  /// Calculates Poynting vector from a moving electron at a point
  ROOT::Math::XYZVector CalcPoyntingVec(const TVector3 fieldPoint, const TVector3 ePosition, const TVector3 eVelocity, const TVector3 eAcceleration);
  
  /// Calculates Poynting vector from electromagnetic field vectors at a point
  /// \param EField Electric field vector
  /// \param BField Magnetic field vector
  /// \Returns The Poynting vector
  ROOT::Math::XYZVector CalcPoyntingVec(const ROOT::Math::XYZVector EField, const ROOT::Math::XYZVector BField);

  /// \brief Calculate the non-relativistic Poynting vector
  /// \param fp Position vector of field point
  /// \param ePos Electron position vector
  /// \param eVel Electron velocity vector (in metres)
  /// \param eAcc Electron acceleration vector
  /// \return Poynting vector in units of W m^-2
  TVector3 CalcPoyntingVecNR(TVector3 fp, TVector3 ePos, TVector3 eVel, TVector3 eAcc);
}

#endif
