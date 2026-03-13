#include "utilities/ConfigParser/SignalConfig.h"

#include "utilities/ConfigParser/YamlUtils.h"

namespace rad {
namespace config {

SignalProcessingConfig ParseSignalProcessingConfig(const YAML::Node& node) {
  SignalProcessingConfig config;
  config.sampleRate =
      GetRequired<double>(node, "sample_rate", "signal_processing");
  config.loFrequency =
      GetRequired<double>(node, "lo_frequency", "signal_processing");

  if (node["acquisition_time"]) {
    config.acquisitionTime = node["acquisition_time"].as<double>();
  }

  if (node["noise"] && node["noise"].IsSequence()) {
    for (size_t i = 0; i < node["noise"].size(); i++) {
      const auto& noiseNode = node["noise"][i];
      NoiseConfig nc;
      std::string ctx = "noise[" + std::to_string(i) + "]";
      nc.temperature = GetRequired<double>(noiseNode, "temperature", ctx);
      nc.resistance = GetRequired<double>(noiseNode, "resistance", ctx);
      nc.seed = GetOptional<int>(noiseNode, "seed", 1234);
      config.noise.push_back(nc);
    }
  }

  return config;
}

OutputConfig ParseOutputConfig(const YAML::Node& node) {
  OutputConfig config;
  config.file = GetRequired<std::string>(node, "file", "output");
  config.group =
      GetOptional<std::string>(node, "group", std::string("/signals"));

  if (node["datasets"] && node["datasets"].IsSequence()) {
    for (const auto& ds : node["datasets"]) {
      config.datasets.push_back(ds.as<std::string>());
    }
  } else {
    config.datasets = {"voltage_i"};
  }

  config.metadata = GetOptional<bool>(node, "metadata", true);
  config.cleanupTrajectories =
      GetOptional<bool>(node, "cleanup_trajectories", true);
  config.trajectoryDir =
      GetOptional<std::string>(node, "trajectory_dir", std::string("/tmp"));

  return config;
}

}  // namespace config
}  // namespace rad
