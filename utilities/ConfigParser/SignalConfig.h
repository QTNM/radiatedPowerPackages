#ifndef SIGNAL_CONFIG_H
#define SIGNAL_CONFIG_H

#include <string>
#include <vector>
#include <yaml-cpp/yaml.h>

namespace rad {
namespace config {

/// Noise source configuration
struct NoiseConfig {
  double temperature;  // Kelvin
  double resistance;   // Ohms
  int seed = 1234;
};

/// Signal processing parameters
struct SignalProcessingConfig {
  double sampleRate;               // Hz
  double loFrequency;              // Hz (not angular)
  std::vector<NoiseConfig> noise;
  double acquisitionTime = -1;     // seconds, -1 means use full trajectory
};

/// Output configuration
struct OutputConfig {
  std::string file;                         // HDF5 output file path
  std::string group = "/signals";           // HDF5 group name
  std::vector<std::string> datasets;        // "voltage_i", "voltage_q", "power_spectrum"
  bool metadata = true;                     // write config metadata as attributes
  bool cleanupTrajectories = true;          // delete temp ROOT files after use
  std::string trajectoryDir = "/tmp";       // directory for temporary trajectory files
};

/// @brief Parse signal processing config from the "signal_processing" YAML node
/// @param node The YAML node
/// @return Populated SignalProcessingConfig
/// @throws std::runtime_error on missing required keys
SignalProcessingConfig ParseSignalProcessingConfig(const YAML::Node& node);

/// @brief Parse output config from the "output" YAML node
/// @param node The YAML node
/// @return Populated OutputConfig
/// @throws std::runtime_error on missing required keys
OutputConfig ParseOutputConfig(const YAML::Node& node);

}  // namespace config
}  // namespace rad

#endif
