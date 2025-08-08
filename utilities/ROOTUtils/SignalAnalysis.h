/// SignalAnalysis.h - ROOT-based signal analysis functions
#ifndef SIGNAL_ANALYSIS_H
#define SIGNAL_ANALYSIS_H

#include "Math/Point3D.h"
#include "Math/Vector3D.h"
#include "TGraph.h"
#include "TVector3.h"

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
                        const ROOT::Math::XYZPoint ePosition, double labTime);

/// @brief Calculate lab time from retarded time (ROOT::Math version)
/// @param fieldPoint Field observation point
/// @param ePosition Electron position
/// @param tRet Retarded time in seconds
/// @return Lab time in seconds
double CalcTimeFromRetardedTime(ROOT::Math::XYZPoint fieldPoint,
                                ROOT::Math::XYZPoint ePosition, double tRet);

/// @brief Calculate lab time from retarded time (TVector3 version)
/// @param fieldPoint Field observation point
/// @param ePosition Electron position
/// @param tRet Retarded time in seconds
/// @return Lab time in seconds
double CalcTimeFromRetardedTime(TVector3 fieldPoint, TVector3 ePosition,
                                double tRet);

}  // namespace rad

#endif