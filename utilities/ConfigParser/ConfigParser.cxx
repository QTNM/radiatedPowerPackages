#include "utilities/ConfigParser/ConfigParser.h"

#include <stdexcept>

#include <yaml-cpp/yaml.h>

#include "utilities/ConfigParser/FieldFactory.h"

namespace rad {
namespace config {

TrajectoryConfig LoadTrajectoryConfig(const std::string& filepath) {
  YAML::Node root;
  try {
    root = YAML::LoadFile(filepath);
  } catch (const YAML::BadFile&) {
    throw std::runtime_error("Cannot open config file: " + filepath);
  } catch (const YAML::ParserException& e) {
    throw std::runtime_error("YAML parse error in " + filepath + ": " +
                             e.what());
  }

  TrajectoryConfig config;

  if (!root["simulation"]) {
    throw std::runtime_error("Config missing required 'simulation' section");
  }
  config.simulation = ParseSimulationConfig(root["simulation"]);

  if (!root["electron"]) {
    throw std::runtime_error("Config missing required 'electron' section");
  }
  config.electron = ParseElectronConfig(root["electron"]);

  if (!root["field"]) {
    throw std::runtime_error("Config missing required 'field' section");
  }
  config.field = CreateField(root["field"]);

  return config;
}

}  // namespace config
}  // namespace rad
