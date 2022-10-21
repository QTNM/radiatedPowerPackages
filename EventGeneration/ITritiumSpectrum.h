/*
  ITritiumSpectrum.h

  Abstract base class for tritium spectrum

  Allows the possiblity for multiple derived classes 
  e.g. are there sterile states? what kind of calculation method are we using?
*/

#ifndef I_TRITIUM_SPECTRUM_H
#define I_TRITIUM_SPECTRUM_H

namespace rad {

  class ITritiumSpectrum {
  public:
    
    /// Function to return the decay rate in s^-1 eV^-1
    /// \param electronEnergy The electron kinetic energy in eV
    /// \return The decay rate in units of s^-1 eV^-1
    virtual double GetDecayRate(double electronEnergy) = 0;

    /// Function for getting random energy according to the distribution
    /// \return An energy in electronvolts
    virtual double DrawRandomEnergy() = 0;
    
  protected:
    double m1; // Mass of first eigenstate
    double m2; // Mass of second eigenstate
    double m3; // Mass of third eigenstate

    // Traditional PMNS mixing angles
    double th12;
    double th13;
    double th23;

    // Spectrum endpoint energy in eV
    double endpoint;

    /// Calculate squared PMNS matrix element Ue1
    /// \return The squared matrix element Ue1
    double CalculateUe1Sq();

    /// Calculate squared PMNS matrix element Ue2
    /// \return The squared matrix element Ue2
    double CalculateUe2Sq();

    /// Calculate squared PMNS matrix element Ue3
    /// \return The squared matrix element Ue3
    double CalculateUe3Sq();
  };

}

#endif
