#ifndef ELECTRON_CONFIG_H
#define ELECTRON_CONFIG_H

#include <cmath>
#include <yaml-cpp/yaml.h>
#include "TVector3.h"

#include "utilities/BasicCore/Constants.h"
#include "utilities/BasicCore/Physics.h"

namespace rad {
namespace config {

/// Electron initial conditions parsed from YAML
struct ElectronConfig {
  TVector3 position;
  TVector3 velocity;
};

/// @brief Compute velocity vector from kinetic energy and pitch angle
/// @param kineticEnergyEV Kinetic energy in eV
/// @param pitchAngleDeg Pitch angle in degrees (angle between velocity and B)
/// @return Velocity 3-vector (perpendicular in x, parallel in z)
inline TVector3 VelocityFromKEAndPitch(double kineticEnergyEV,
                                       double pitchAngleDeg) {
  double pitchRad = pitchAngleDeg * PI / 180.0;
  double speed = GetSpeedFromKE(kineticEnergyEV, ME);
  return TVector3(speed * std::sin(pitchRad), 0.0, speed * std::cos(pitchRad));
}

/// @brief Parse electron config from a YAML node
///
/// Supports two mutually exclusive modes:
///   - kinetic_energy + pitch_angle: computes velocity relativistically
///   - velocity: uses the raw 3-vector directly
///
/// @param node The "electron" YAML node
/// @return Populated ElectronConfig with position and velocity
/// @throws std::runtime_error on missing/conflicting fields
ElectronConfig ParseElectronConfig(const YAML::Node& node);

}  // namespace config
}  // namespace rad

#endif
