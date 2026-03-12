#ifndef SIMULATION_CONFIG_H
#define SIMULATION_CONFIG_H

#include <limits>
#include <string>
#include <yaml-cpp/yaml.h>

namespace rad {
namespace config {

/// Geometric bounds for early trajectory termination
struct GeometricBounds {
  double zMin = -std::numeric_limits<double>::max();
  double zMax = std::numeric_limits<double>::max();
  double rMax = std::numeric_limits<double>::max();
};

/// Simulation parameters parsed from YAML
struct SimulationConfig {
  std::string outputFile;
  double time;
  double stepSize;
  bool energyLoss = true;
  double initialTime = 0.0;
  GeometricBounds bounds;
};

/// @brief Parse simulation config from a YAML node
/// @param node The "simulation" YAML node
/// @return Populated SimulationConfig
/// @throws std::runtime_error on missing required fields
SimulationConfig ParseSimulationConfig(const YAML::Node& node);

}  // namespace config
}  // namespace rad

#endif
