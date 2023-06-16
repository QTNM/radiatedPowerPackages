/*
  HalfWaveDipole.h
  Hypothetical antenna with a length of half a wavelength
*/

#ifndef HALF_WAVE_DIPOLE_H
#define HALF_WAVE_DIPOLE_H

#include "Antennas/IAntenna.h"

#include "TVector3.h"

namespace rad {

  class HalfWaveDipole : public IAntenna {
    
  public:

    /// \param antPos The position of the antenna
    /// \param antXAx Antenna defined X axis
    /// \param antZAx Antenna defined Z axis
    /// \param freq Central frequency of the antenna
    /// \param delay Time delay of this antenna (in seconds)
    HalfWaveDipole(TVector3 antPos, TVector3 antXAx, TVector3 antZAx, double freq, double delay=0.0);
    
    TVector3 GetETheta(const TVector3 electronPosition);
    TVector3 GetEPhi(const TVector3 electronPosition);

    double GetETheta(double theta, double phi);

    double GetEPhi(double theta, double phi);

    double GetHEff();

    double GetHEff(TVector3 ePos);

    /// Calculates the effective area (position dependent) of the dipole
    /// \param ePos Radiating electron position vector
    /// \return The effective area in metres squared
    double GetAEff(TVector3 ePos); 

    double GetAEffTheta(TVector3 ePos);

    double GetAEffPhi(TVector3 ePos);

  private:
    double PRad; // Surface integral of radiation pattern, used for normalisation

  };
}

#endif
