/*
  IsotropicAntenna.h

  A model of an antenna with a completely uniform radiation pattern.
  The effective area of this antenna can be arbitrarily scaled
*/

#ifndef ISOTROPIC_ANTENNA_H
#define ISOTROPIC_ANTENNA_H

#include "Antennas/IAntenna.h"

namespace rad {
class IsotropicAntenna : public IAntenna {
 public:
  /// @brief Constructor
  /// @param antPos Antenna position vector
  /// @param AEff Chosen effective area in m^2
  /// @param Z Chosen antenna impedance in Ohms
  /// @param freq
  /// @param delay
  IsotropicAntenna(TVector3 antPos, double AEff, double Z, double freq = 27e9,
                   double delay = 0);

  //////////// Radiation pattern functions ////////////

  /// @brief Electric field vector in theta direction
  /// @param ePos Electron position vector
  /// @return Electric field vector
  TVector3 GetETheta(const TVector3 ePos) override;

  /// @brief Electric field vector in phi direction
  /// @param ePos Electron position vector
  /// @return Electric field vector
  TVector3 GetEPhi(const TVector3 ePos) override;

  /// @brief Electric field vector in theta direction
  /// @param theta Theta angle
  /// @param phi Phi angle
  /// @return Electric field vector
  double GetETheta(double theta, double phi) override;

  /// @brief Electric field vector in phi direction
  /// @param theta Theta angle
  /// @param phi Phi angle
  /// @return Electric field vector
  double GetEPhi(double theta, double phi) override;

  /// @brief Gets the effective antenna height/length
  /// @return Effective height in metres
  double GetHEff() override;

  /// @brief Gets the effective antenna height/length
  /// @param ePos Electron position vector
  /// @return Effective height in metres
  double GetHEff(TVector3 ePos) override;

  /// @brief Virtual function for returning antenna effective area
  /// @param electronPosition The position vector of the electron
  /// @return The effective area of the antenna (in metres)
  double GetAEff(TVector3 electronPosition) override;

  double GetAEffTheta(TVector3 electronPosition) override;

  double GetAEffPhi(TVector3 electronPosition) override;

 private:
  double effectiveArea;
  double imp;
  double PRad;
};
}  // namespace rad

#endif