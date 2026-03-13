#ifndef PIPELINE_CONFIG_H
#define PIPELINE_CONFIG_H

#include <memory>
#include <string>
#include <vector>

#include "physics/ElectronDynamics/BaseField.h"
#include "utilities/ConfigParser/DetectorFactory.h"
#include "utilities/ConfigParser/ElectronConfig.h"
#include "utilities/ConfigParser/SignalConfig.h"
#include "utilities/ConfigParser/SimulationConfig.h"

namespace rad {
namespace config {

/// Full pipeline configuration parsed from a YAML file
struct PipelineConfig {
  SimulationConfig simulation;
  std::vector<ElectronConfig> electrons;
  std::unique_ptr<BaseField> field;
  Detector detector;
  SignalProcessingConfig signalProcessing;
  OutputConfig output;
};

/// @brief Load and parse a full pipeline configuration file
/// @param filepath Path to the YAML configuration file
/// @return Fully populated PipelineConfig
/// @throws std::runtime_error on parse errors or missing sections
PipelineConfig LoadPipelineConfig(const std::string& filepath);

}  // namespace config
}  // namespace rad

#endif
