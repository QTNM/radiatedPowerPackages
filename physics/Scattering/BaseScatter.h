/*
  BaseScatter.h

  Abstract base class for calculating scattering cross-sections
*/

#ifndef BASE_SCATTER_H
#define BASE_SCATTER_H

#include "TVector3.h"

namespace rad {
class BaseScatter {
 private:
  double ke;  // Kinetic energy of particle in eV

 protected:
  /// @brief Getter function for incident electron KE
  /// @return The incident electron KE in eV
  double GetIncidentKE() { return ke; }

 public:
  /// @brief Parametrised constructor
  /// @param T Incident electron ke in eV
  BaseScatter(double T) : ke(T) {}

  /// @brief Virtual function for getting the total xsec
  /// @return Total cross-section in m^2
  virtual double GetTotalXSec() = 0;

  /// @brief Calculates the mean free path for a given xsec
  /// @param N Number density of target atoms
  /// @return Mean free path in metres
  double GetMeanFreePath(double N);

  /// @brief Calculates particle mean free time before scattering
  /// @param N Number density of atoms/molecules
  /// @return Mean free time in seconds
  double GetMeanFreeTime(double N);

  /// @brief Calculate random vector after scattering
  /// @param vel Unscattered velocity vector
  /// @param outKE Kinetic energy of outgoing particle in eV
  /// @param theta Known scattering angle in radians
  /// @return Scattered velocity vector
  TVector3 GetScatteredVector(TVector3 vel, double outKE, double theta);

  /// @brief Sets incident kinetic energy
  /// @param T New incident KE in eV
  void SetIncidentKE(double T) { ke = T; }

  /// @brief Virtual destructor
  virtual ~BaseScatter(){};
};
}  // namespace rad

#endif