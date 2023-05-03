/*
  BaseScatter.h

  Abstract base class for calculating scattering cross-sections
*/

#ifndef BASE_SCATTER_H
#define BASE_SCATTER_H

namespace rad {
class BaseScatter {
 private:
  double ke;  // Kinetic energy of particle in eV

 protected:
  /// @brief Getter function for incident electron KE
  /// @return
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

  /// @brief
  /// @param N
  /// @return
  double GetMeanFreeTime(double N);

  virtual ~BaseScatter(){};
};
}  // namespace rad

#endif