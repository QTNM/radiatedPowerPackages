/*
  HertzianDipole.h
  Simplest case of a Hertzian dipole
*/

#ifndef HERTZIAN_DIPOLE_H
#define HERTZIAN_DIPOLE_H

#include "Antennas/IAntenna.h"

#include "TVector3.h"

namespace rad {

  class HertzianDipole : public IAntenna {
    
  public:
    
    /// \param antPos The position of the antenna
    /// \param antXAx Antenna defined X axis
    /// \param antZAx Antenna defined Z axis
    /// \param freq Central frequency of the antenna
    /// \param delay Time delay of this antenna
    HertzianDipole(TVector3 antPos, TVector3 antXAx, TVector3 antZAx, double freq, double delay=0.0);
    
    TVector3 GetETheta(const TVector3 electronPosition);
    TVector3 GetEPhi(const TVector3 electronPosition);

    double GetETheta(double theta, double phi);

    double GetEPhi(double theta, double phi);

    double GetHEff();

    double GetHEff(TVector3 ePos);

    /// Returns the effective area of a Hertzian dipole
    /// \param ePos Electron position vector (in metres)
    /// \return The antenna effective area (in metres squared)
    double GetAEff(TVector3 ePos);

    double GetAEffTheta(TVector3 ePos);

    double GetAEffPhi(TVector3 ePos);
  };
}

#endif
