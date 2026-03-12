#ifndef ELECTRON_CONFIG_H
#define ELECTRON_CONFIG_H

#include <yaml-cpp/yaml.h>
#include "TVector3.h"

namespace rad {
namespace config {

/// Electron initial conditions parsed from YAML
struct ElectronConfig {
  TVector3 position;
  TVector3 velocity;
};

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
