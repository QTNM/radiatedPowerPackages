/*
  BasicTritiumSpectrum.cxx

  Decay electron energy spectrum calculated in a rudimentary fashion
  See Eur. Phys. J. C (2019) 79:204 for further details
*/

#ifndef BASIC_TRITIUM_SPECTRUM_H
#define BASIC_TRITIUM_SPECTRUM_H

#include "EventGeneration/ITritiumSpectrum.h"

#include "TF1.h"

namespace rad
{
  class BasicTritiumSpectrum : public ITritiumSpectrum {
  public:
    /// Constructor for the spectrum
    /// \param mass1 Mass of the first eigenstate in eV c^-2
    /// \param mass2 Mass of the second eigenstate in eV c^-2
    /// \param mass3 Mass of the third eigenstate in eV c^-2
    /// \param theta12 Value of theta12 in radians
    /// \param theta13 Value of theta13 in radians
    /// \param theta23 Value of theta23 in radians
    /// \param endpointE Value of the endpoint in eV
    BasicTritiumSpectrum(double mass1, double mass2, double mass3, double theta12, double theta13, double theta23, double endpointE);

    /// Calculates the decay rate of a tritium atom for a given decay electron kinetic energy
    /// \param electronEnergy The kinetic energy of the decay electron (in eV)
    /// \Returns The decay rate with units of s^-1 eV^-1
    double GetDecayRate(double electronEnergy);

    /// Draws a random kinetic energy from the distribution
    /// \Returns a random energy between 0 and the endpoint energy in eV
    double DrawRandomEnergy();

  private:
    double Ue1Sq;
    double Ue2Sq;
    double Ue3Sq;
    double maxRate;
  };
}

#endif
