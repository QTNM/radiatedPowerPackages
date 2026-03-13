#ifndef SCATTERING_CONFIG_H
#define SCATTERING_CONFIG_H

#include <string>
#include <vector>
#include <yaml-cpp/yaml.h>

namespace rad {
namespace config {

/// A single gas species with its number density
struct GasSpecies {
  std::string species;  // "H", "H2", "He", "T", "T2"
  double density;       // number density in m^-3
};

/// Scattering configuration: one or more gas species
struct ScatteringConfig {
  std::vector<GasSpecies> gases;
};

/// @brief Parse scattering config from the "scattering" YAML node
///
/// Expects a "gases" sequence where each entry has:
///   species: H | H2 | He | T | T2
///   density: <double>   (m^-3)
///
/// @param node The "scattering" YAML node
/// @return Populated ScatteringConfig
/// @throws std::runtime_error on missing required fields or unknown species
ScatteringConfig ParseScatteringConfig(const YAML::Node& node);

}  // namespace config
}  // namespace rad

#endif
