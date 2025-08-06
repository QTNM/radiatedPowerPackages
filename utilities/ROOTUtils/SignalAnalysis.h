/// SignalAnalysis.h - ROOT-based signal analysis functions
#ifndef SIGNAL_ANALYSIS_H 
#define SIGNAL_ANALYSIS_H

#include "TGraph.h"
#include "TVector3.h"
#include "Math/Point3D.h"
#include "Math/Vector3D.h"

namespace rad {

/// @brief Calculate effective area of Hertzian dipole antenna
/// @param wavelength Signal wavelength in meters
/// @param dipoleDir Dipole direction vector
/// @param ePosition Electron position
/// @param antennaPoint Antenna position  
/// @return Effective area in square meters
double CalcAeHertzianDipole(double wavelength,
                            const ROOT::Math::XYZVector dipoleDir,
                            const ROOT::Math::XYZPoint ePosition,
                            const ROOT::Math::XYZPoint antennaPoint);

/// @brief Calculate effective length of Hertzian dipole antenna
/// @param wavelength Signal wavelength in meters
/// @param dipoleDir Dipole direction vector
/// @param ePosition Electron position
/// @param antennaPoint Antenna position
/// @return Effective length in meters
double CalcAlHertzianDipole(double wavelength,
                            const ROOT::Math::XYZVector dipoleDir,
                            const ROOT::Math::XYZPoint ePosition,
                            const ROOT::Math::XYZPoint antennaPoint);

/// @brief Calculate retarded time for electromagnetic propagation
/// @param fieldPoint Field observation point
/// @param ePosition Electron position at emission
/// @param labTime Lab time of emission
/// @return Retarded time in seconds
double CalcRetardedTime(const ROOT::Math::XYZPoint fieldPoint,
                        const ROOT::Math::XYZPoint ePosition,
                        double labTime);

/// @brief Calculate lab time from retarded time (ROOT::Math version)
/// @param fieldPoint Field observation point
/// @param ePosition Electron position
/// @param tRet Retarded time in seconds
/// @return Lab time in seconds
double CalcTimeFromRetardedTime(ROOT::Math::XYZPoint fieldPoint,
                                ROOT::Math::XYZPoint ePosition, 
                                double tRet);

/// @brief Calculate lab time from retarded time (TVector3 version)
/// @param fieldPoint Field observation point
/// @param ePosition Electron position  
/// @param tRet Retarded time in seconds
/// @return Lab time in seconds
double CalcTimeFromRetardedTime(TVector3 fieldPoint, TVector3 ePosition, double tRet);

/// @brief Get particle speed from kinetic energy
/// @param T Particle kinetic energy in eV
/// @param particleMass Particle mass in kg
/// @return Speed in m/s
double GetSpeedFromKE(double T, double particleMass);

/// @brief Calculate gyroradius/Larmor radius of charged particle
/// @param velocity Velocity vector in m/s
/// @param bField Magnetic field vector in Tesla
/// @param particleMass Particle mass in kg
/// @return Gyroradius in meters
double GetGyroradius(TVector3 velocity, TVector3 bField, double particleMass);

/// @brief Calculate relativistic cyclotron frequency
/// @param BField Magnetic field vector in Tesla
/// @param charge Particle charge in Coulombs (default: electron)
/// @param energy Kinetic energy in eV
/// @param mass Particle mass in kg (default: electron)
/// @return Angular frequency vector in rad/s
TVector3 calculate_omega(const TVector3 BField,
                         double charge = -1.602176634e-19,  // electron charge
                         double energy = 0.0,
                         double mass = 9.1093837015e-31);   // electron mass

/// @brief Calculate relativistic electron cyclotron frequency
/// @param KE Electron kinetic energy in eV  
/// @param B Magnetic field strength in Tesla
/// @return Cyclotron frequency in Hz
double CalcCyclotronFreq(double KE, double B = 1.0);

}  // namespace rad

#endif