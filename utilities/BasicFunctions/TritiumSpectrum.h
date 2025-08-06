/*
  TritiumSpectrum.h

  Contains functions relating to the tritium spectra
*/
#ifndef TRITIUM_SPECTRUM_H
#define TRITIUM_SPECTRUM_H

namespace rad {
/// Gives the differential decay rate of a tritium nucleus
/// See Eur. Phys. J. C (2019) 79:204 for further details
/// \param electronT Electron kinetic energy (in eV)
/// \param m1 Mass of the m1 eigenstate (in eV/c^2)
/// \param m2 Mass of the m2 eigenstate (in eV/c^2)
/// \param m3 Mass of the m3 eigenstate (in eV/c^2)
/// \param endpointE The inputted tritium endpoint energy (in eV)
/// Returns the decay rate in s^-1 eV^-1
double TritiumDecayRate(double electronT, double m1, double m2, double m3,
                        double endpointE = 18574);
}  // namespace rad

#endif
