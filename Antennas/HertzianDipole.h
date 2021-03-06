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

    double GetHEff();
  };
}

#endif
