#include "utilities/ConfigParser/ElectronConfig.h"

#include <stdexcept>

#include "utilities/ConfigParser/YamlUtils.h"

namespace rad {
namespace config {

ElectronConfig ParseElectronConfig(const YAML::Node& node) {
  ElectronConfig config;

  if (!node["position"]) {
    throw std::runtime_error("electron: missing required key 'position'");
  }
  config.position = ParseVector3(node["position"], "electron.position");

  bool hasKE = node["kinetic_energy"].IsDefined();
  bool hasVel = node["velocity"].IsDefined();

  if (hasKE && hasVel) {
    throw std::runtime_error(
        "electron: specify either kinetic_energy+pitch_angle or velocity, "
        "not both");
  }
  if (!hasKE && !hasVel) {
    throw std::runtime_error(
        "electron: must specify either kinetic_energy+pitch_angle or velocity");
  }

  if (hasKE) {
    double ke = node["kinetic_energy"].as<double>();
    if (!node["pitch_angle"]) {
      throw std::runtime_error(
          "electron: kinetic_energy requires pitch_angle");
    }
    double pitchDeg = node["pitch_angle"].as<double>();
    config.velocity = VelocityFromKEAndPitch(ke, pitchDeg);
  } else {
    config.velocity = ParseVector3(node["velocity"], "electron.velocity");
  }

  return config;
}

}  // namespace config
}  // namespace rad
