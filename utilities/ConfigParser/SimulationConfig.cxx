#include "utilities/ConfigParser/SimulationConfig.h"

#include <stdexcept>

namespace rad {
namespace config {

SimulationConfig ParseSimulationConfig(const YAML::Node& node) {
  SimulationConfig config;

  if (node["output_file"]) {
    config.outputFile = node["output_file"].as<std::string>();
  }

  if (!node["time"]) {
    throw std::runtime_error("simulation: missing required key 'time'");
  }
  config.time = node["time"].as<double>();

  if (!node["step_size"]) {
    throw std::runtime_error("simulation: missing required key 'step_size'");
  }
  config.stepSize = node["step_size"].as<double>();

  if (node["energy_loss"]) {
    config.energyLoss = node["energy_loss"].as<bool>();
  }

  if (node["initial_time"]) {
    config.initialTime = node["initial_time"].as<double>();
  }

  // Parse optional geometric bounds
  if (node["bounds"]) {
    const YAML::Node& b = node["bounds"];
    if (b["z_min"]) {
      config.bounds.zMin = b["z_min"].as<double>();
    }
    if (b["z_max"]) {
      config.bounds.zMax = b["z_max"].as<double>();
    }
    if (b["r_max"]) {
      config.bounds.rMax = b["r_max"].as<double>();
    }
  }

  return config;
}

}  // namespace config
}  // namespace rad
