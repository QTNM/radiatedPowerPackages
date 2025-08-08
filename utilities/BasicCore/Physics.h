// Physics.h
// Header-only inline functions for basic physics calculations

#ifndef PHYSICS_H
#define PHYSICS_H

#include <cmath>

#include "Constants.h"

namespace rad {
/// @brief Calculate retarded time for electromagnetic field calculations
/// @param fieldPointPos Field point position (x, y, z) in meters
/// @param ePositionPos Source position (x, y, z) in meters
/// @param labTime Lab time in seconds
/// @return Retarded time in seconds
inline double CalcRetardedTime(const double fieldPointPos[3],
                               const double ePositionPos[3],
                               const double labTime) {
  double dx = ePositionPos[0] - fieldPointPos[0];
  double dy = ePositionPos[1] - fieldPointPos[1];
  double dz = ePositionPos[2] - fieldPointPos[2];
  double distance = std::sqrt(dx * dx + dy * dy + dz * dz);
  return labTime - distance / C;
}

/// @brief Calculate lab time from retarded time
/// @param fieldPointPos Field point position (x, y, z) in meters
/// @param ePositionPos Source position (x, y, z) in meters
/// @param tRet Retarded time in seconds
/// @return Lab time in seconds
inline double CalcTimeFromRetardedTime(const double fieldPointPos[3],
                                       const double ePositionPos[3],
                                       double tRet) {
  double dx = ePositionPos[0] - fieldPointPos[0];
  double dy = ePositionPos[1] - fieldPointPos[1];
  double dz = ePositionPos[2] - fieldPointPos[2];
  double distance = std::sqrt(dx * dx + dy * dy + dz * dz);
  return tRet + distance / C;
}

/// @brief Get particle speed from kinetic energy
/// @param T Particle kinetic energy in eV
/// @param particleMass Particle mass in kg
/// @return Speed in m/s
inline double GetSpeedFromKE(double T, double particleMass) {
  double gamma = T * QE / (particleMass * C * C) + 1;
  double betaSq = 1 - 1 / (gamma * gamma);
  double speed = sqrt(betaSq) * C;
  return speed;
}

/// @brief Calculate the free space radiated power
/// @param ke Kinetic energy in eV
/// @param B Magnetic field strength in tesla
/// @param theta Pitch angle in radians
/// @param m Particle mass in kg
/// @return Larmor power in Watts
inline double CalcLarmorPower(double ke, double B, double theta,
                              double m = ME) {
  const double f0{QE * B / (m * (2 * PI))};
  const double beta{GetSpeedFromKE(ke, m) / C};
  return (2 * PI) * pow(QE * f0 * beta * sin(theta), 2) /
         ((3 * EPSILON0 * C) * (1 - beta * beta));
}

/// @brief Function for a chirp signal
/// @param A Signal amplitude
/// @param t Time at which to generate signal [s]
/// @param phi0 Initial phase [radians]
/// @param f0 Initial frequency [Hz]
/// @param c Chirp rate [Hz s^-1]
/// @return The chirp function at the supplied time
inline double ChirpSignal(double A, double t, double phi0, double f0,
                          double c) {
  return A * sin(phi0 + 2 * PI * (c * t * t / 2 + f0 * t));
}

/// @brief Calculate relativistic electron cyclotron frequency
/// @param KE Electron kinetic energy in eV
/// @param B Magnetic field strength in Tesla
/// @return Cyclotron frequency in Hz
inline double CalcCyclotronFreq(double KE, double B = 1.0) {
  double freq{QE * B / (ME + (KE * QE / pow(C, 2)))};
  return freq / (2 * PI);
}

}  // namespace rad

#endif