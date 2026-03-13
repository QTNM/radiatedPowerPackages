#ifndef CONFIG_PARSER_H
#define CONFIG_PARSER_H

#include <memory>
#include <string>
#include <vector>

#include "physics/ElectronDynamics/BaseField.h"
#include "utilities/ConfigParser/ElectronConfig.h"
#include "utilities/ConfigParser/SimulationConfig.h"

namespace rad {
namespace config {

/// Full trajectory configuration parsed from a YAML file
struct TrajectoryConfig {
  SimulationConfig simulation;
  std::vector<ElectronConfig> electrons;
  std::unique_ptr<BaseField> field;
};

/// @brief Load and parse a trajectory configuration file
/// @param filepath Path to the YAML configuration file
/// @return Fully populated TrajectoryConfig
/// @throws std::runtime_error on parse errors or missing sections
TrajectoryConfig LoadTrajectoryConfig(const std::string& filepath);

}  // namespace config
}  // namespace rad

#endif
