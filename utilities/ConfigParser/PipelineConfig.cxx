#include "utilities/ConfigParser/PipelineConfig.h"

#include <stdexcept>

#include <yaml-cpp/yaml.h>

#include "utilities/ConfigParser/ElectronGenerator.h"
#include "utilities/ConfigParser/FieldFactory.h"

namespace rad {
namespace config {

PipelineConfig LoadPipelineConfig(const std::string& filepath) {
  YAML::Node root;
  try {
    root = YAML::LoadFile(filepath);
  } catch (const YAML::BadFile&) {
    throw std::runtime_error("Cannot open config file: " + filepath);
  } catch (const YAML::ParserException& e) {
    throw std::runtime_error("YAML parse error in " + filepath + ": " +
                             e.what());
  }

  PipelineConfig config;

  // Required: simulation
  if (!root["simulation"]) {
    throw std::runtime_error("Config missing required 'simulation' section");
  }
  config.simulation = ParseSimulationConfig(root["simulation"]);

  // Required: electron(s)
  bool hasSingular = root["electron"].IsDefined();
  bool hasPlural = root["electrons"].IsDefined();

  if (hasSingular && hasPlural) {
    throw std::runtime_error(
        "Config has both 'electron' and 'electrons' — use one or the other");
  }
  if (!hasSingular && !hasPlural) {
    throw std::runtime_error(
        "Config missing required 'electron' or 'electrons' section");
  }

  if (hasSingular) {
    config.electrons.push_back(ParseElectronConfig(root["electron"]));
  } else {
    config.electrons = GenerateElectrons(root["electrons"]);
  }

  // Required: field
  if (!root["field"]) {
    throw std::runtime_error("Config missing required 'field' section");
  }
  config.field = CreateField(root["field"]);

  // Optional: scattering
  if (root["scattering"]) {
    config.scattering = ParseScatteringConfig(root["scattering"]);
  }

  // Required: detector
  if (!root["detector"]) {
    throw std::runtime_error("Config missing required 'detector' section");
  }
  config.detector = CreateDetector(root["detector"]);

  // Required: signal_processing
  if (!root["signal_processing"]) {
    throw std::runtime_error(
        "Config missing required 'signal_processing' section");
  }
  config.signalProcessing =
      ParseSignalProcessingConfig(root["signal_processing"]);

  // Required: output
  if (!root["output"]) {
    throw std::runtime_error("Config missing required 'output' section");
  }
  config.output = ParseOutputConfig(root["output"]);

  return config;
}

}  // namespace config
}  // namespace rad
