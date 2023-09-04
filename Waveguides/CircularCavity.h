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
  /// @param probe Desired probe position
  CircularCavity(double radius, double length, TVector3 probe);

  /// @brief Getter function for radius
  /// @return Cavity inner radius in metres
  double GetRadius() { return a; }

  /// @brief Getter function for length
  /// @return Cavity length in metres
  double GetLength() { return d; }

  double GetVolume() { return TMath::Pi() * a * a * d; }

  /// @brief Calculate resonant frequency of given mode
  /// @param modeType TE or TM
  /// @param n Angular mode number
  /// @param m Radial mode number
  /// @param l Axial mode index
  /// @return Frequency in Hertz
  double GetResonantModeF(Mode_t modeType, unsigned int n, unsigned int m,
                          unsigned int l) override;

  /// @brief Gets the cutoff frequency for a particular waveguide mode
  /// @param modeType The mode type to get (TE, TM, TEM)
  /// @param n The angular mode number
  /// @param m The radial mode number
  /// @return The cutoff frequency in Hertz
  double GetCutoffFrequency(Mode_t modeType, int n, int m) override;

  /// @brief Calculate electric field for given mode
  /// @param rho Radial coordinate in metres
  /// @param phi Angular coordinate in radians
  /// @param z Axial coordinate in metres
  /// @param modeType TE or TM
  /// @param A Constant factor
  /// @param m Angular mode number
  /// @param n Radial mode number
  /// @param p Axial mode number
  /// @param t Time in seconds
  /// @param state Choose polarisation state
  /// @return Complex 3-vector of E field
  ComplexVector3 GetModeEField(double rho, double phi, double z,
                               Mode_t modeType, double A, unsigned int m,
                               unsigned int n, unsigned int p,
                               bool state = true, double t = 0);

  /// @brief Calculate electric field for given mode
  /// @param pos Position vector (in metres)
  /// @param modeType TE or TM
  /// @param A Constant factor
  /// @param m Angular mode number
  /// @param n Radial mode number
  /// @param p Axial mode number
  /// @param t Time in seconds
  /// @param state Choose polarisation state
  /// @return Complex 3-vector of E field
  ComplexVector3 GetModeEField(TVector3 pos, Mode_t modeType, double A,
                               unsigned int m, unsigned int n, unsigned int p,
                               bool state = true, double t = 0);

  /// @brief Get mode electric field with no time variation
  /// @param pos Position to evaluate field at
  /// @param modeType Select TE or TM
  /// @param A Normalisation factor
  /// @param m Angular mode number
  /// @param n Radial mode number
  /// @param p Axial mode number
  /// @param state Choose from two polarisation states
  /// @return Complex vector of electric field
  ComplexVector3 GetModalEField(TVector3 pos, Mode_t modeType, double A,
                                unsigned int m, unsigned int n, unsigned int p,
                                bool state) override;

  /// @brief Calculate electric field at maximum time
  /// @param pos Position 3-vector (in metres)
  /// @param modeType TE or TM
  /// @param A Constant factor
  /// @param m Angular mode number
  /// @param n Radial mode number
  /// @param p Axial mode number
  /// @param state Choose polarisation state
  /// @return Complex 3-vector or E field
  ComplexVector3 GetMaxEField(TVector3 pos, Mode_t modeType, double A,
                              unsigned int m, unsigned int n, unsigned int p,
                              bool state);

  /// @brief Calculate H field for given mode
  /// @param rho Radial coordinate in metres
  /// @param phi Angular coordinate in radians
  /// @param z Axial coordinate in metres
  /// @param modeType TE or TM
  /// @param A Constant factor
  /// @param m Angular mode number
  /// @param n Radial mode number
  /// @param p Axial mode number
  /// @param state Polarisation state
  /// @param t Time in seconds
  /// @return Complex 3-vector of H field
  ComplexVector3 GetModeHField(double rho, double phi, double z,
                               Mode_t modeType, double A, unsigned int m,
                               unsigned int n, unsigned int p,
                               bool state = true, double t = 0);

  /// @brief Calculate H field for given mode
  /// @param pos Position vector (in metres)
  /// @param modeType TE or TM
  /// @param A Constant factor
  /// @param m Angular mode number
  /// @param n Radial mode number
  /// @param p Axial mode number
  /// @param state Polarisation state
  /// @param t Time in seconds
  /// @return Complex 3-vector of H field
  ComplexVector3 GetModeHField(TVector3 pos, Mode_t modeType, double A,
                               unsigned int m, unsigned int n, unsigned int p,
                               bool state = true, double t = 0);
};
}  // namespace rad

#endif