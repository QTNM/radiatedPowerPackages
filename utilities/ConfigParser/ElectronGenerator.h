#ifndef ELECTRON_GENERATOR_H
#define ELECTRON_GENERATOR_H

#include <vector>
#include <yaml-cpp/yaml.h>

#include "utilities/ConfigParser/ElectronConfig.h"

namespace rad {
namespace config {

/// @brief Generate multiple electron initial conditions from a YAML node
///
/// kinetic_energy can be:
///   - A fixed scalar value (same for all electrons)
///   - A scan: {min, max} (linearly spaced over count, inclusive)
///   - A uniform: {min, max} (uniformly distributed random values)
///
/// Velocity direction is set by exactly one of:
///   - pitch_angle: fixed, scan, or uniform (velocity in x-z plane)
///   - isotropic: true (uniform random direction on the unit sphere)
///
/// Position can be:
///   - A fixed [x, y, z] array
///   - disk: {radius, z} (uniform random on a disk in x-y at given z)
///   - cylinder: {radius, z_min, z_max} (uniform random in a cylinder)
///
/// @param node The "electrons" YAML node (must contain "count")
/// @return Vector of ElectronConfig, one per electron
/// @throws std::runtime_error on invalid configuration
std::vector<ElectronConfig> GenerateElectrons(const YAML::Node& node);

}  // namespace config
}  // namespace rad

#endif
