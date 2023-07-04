/*
  Circular resonant cavity

  S. Jones 30-06-2023
*/

#ifndef CIRCULAR_CAVITY_H
#define CIRCULAR_CAVITY_H

#include "TMath.h"
#include "Waveguides/ICavity.h"

namespace rad {
class CircularCavity : public ICavity {
 private:
  double a;
  double d;

 public:
  /// @brief Parametrised constructor
  /// @param radius Cavity radius in metres
  /// @param length Cavity length in metres
  CircularCavity(double radius, double length) : a(radius), d(length) {}

  /// @brief Calculate resonant frequency of given mode
  /// @param modeType TE or TM
  /// @param n Angular mode number
  /// @param m Radial mode number
  /// @param l Axial mode index
  /// @return Frequency in Hertz
  double GetResonantModeF(Mode_t modeType, unsigned int n, unsigned int m,
                          unsigned int l) override;

  /// @brief Calculate electric field for given mode
  /// @param rho Radial coordinate in metres
  /// @param phi Angular coordinate in radians
  /// @param z Axial coordinate in metres
  /// @param modeType TE or TM
  /// @param A Constant factor
  /// @param n Angular mode number
  /// @param m Radial mode number
  /// @param l Axial mode number
  /// @return Complex 3-vector of E field
  ComplexVector3 GetModeEField(double rho, double phi, double z,
                               Mode_t modeType, double A, unsigned int n,
                               unsigned int m, unsigned int l);

  /// @brief Calculate H field for given mode
  /// @param rho Radial coordinate in metres
  /// @param phi Angular coordinate in radians
  /// @param z Axial coordinate in metres
  /// @param modeType TE or TM
  /// @param A Constant factor
  /// @param n Angular mode number
  /// @param m Radial mode number
  /// @param l Axial mode number
  /// @return Complex 3-vector of H field
  ComplexVector3 GetModeHField(double rho, double phi, double z,
                               Mode_t modeType, double A, unsigned int n,
                               unsigned int m, unsigned int l);
};
}  // namespace rad

#endif