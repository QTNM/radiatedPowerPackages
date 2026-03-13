#include "utilities/ConfigParser/ElectronGenerator.h"

#include <cmath>
#include <random>
#include <stdexcept>
#include <string>

#include "utilities/BasicCore/Constants.h"
#include "utilities/BasicCore/Physics.h"
#include "utilities/ConfigParser/YamlUtils.h"

namespace rad {
namespace config {

/// Resolve N values for a scalar parameter (kinetic_energy or pitch_angle)
static std::vector<double> ResolveScalarParam(const YAML::Node& node,
                                              const std::string& name,
                                              int count,
                                              std::mt19937& rng) {
  std::vector<double> values(count);

  if (node.IsScalar()) {
    double val = node.as<double>();
    std::fill(values.begin(), values.end(), val);
  } else if (node.IsMap()) {
    if (node["scan"]) {
      double min = GetRequired<double>(node["scan"], "min", name + ".scan");
      double max = GetRequired<double>(node["scan"], "max", name + ".scan");
      for (int i = 0; i < count; i++) {
        values[i] =
            count > 1 ? min + (max - min) * double(i) / double(count - 1)
                      : min;
      }
    } else if (node["uniform"]) {
      double min =
          GetRequired<double>(node["uniform"], "min", name + ".uniform");
      double max =
          GetRequired<double>(node["uniform"], "max", name + ".uniform");
      std::uniform_real_distribution<double> dist(min, max);
      for (int i = 0; i < count; i++) {
        values[i] = dist(rng);
      }
    } else {
      throw std::runtime_error("electrons." + name +
                               ": expected a scalar, scan, or uniform");
    }
  } else {
    throw std::runtime_error("electrons." + name +
                             ": expected a scalar or map");
  }

  return values;
}

/// Resolve N positions
static std::vector<TVector3> ResolvePositions(const YAML::Node& node,
                                              int count, std::mt19937& rng) {
  std::vector<TVector3> positions(count);

  if (node.IsSequence()) {
    TVector3 fixed = ParseVector3(node, "electrons.position");
    std::fill(positions.begin(), positions.end(), fixed);
  } else if (node.IsMap()) {
    if (node["disk"]) {
      double radius =
          GetRequired<double>(node["disk"], "radius", "position.disk");
      double z = GetRequired<double>(node["disk"], "z", "position.disk");
      std::uniform_real_distribution<double> uDist(0.0, 1.0);
      std::uniform_real_distribution<double> angleDist(0.0, 2.0 * PI);
      for (int i = 0; i < count; i++) {
        double r = radius * std::sqrt(uDist(rng));
        double theta = angleDist(rng);
        positions[i] = TVector3(r * std::cos(theta), r * std::sin(theta), z);
      }
    } else if (node["cylinder"]) {
      double radius =
          GetRequired<double>(node["cylinder"], "radius", "position.cylinder");
      double zMin =
          GetRequired<double>(node["cylinder"], "z_min", "position.cylinder");
      double zMax =
          GetRequired<double>(node["cylinder"], "z_max", "position.cylinder");
      std::uniform_real_distribution<double> uDist(0.0, 1.0);
      std::uniform_real_distribution<double> angleDist(0.0, 2.0 * PI);
      std::uniform_real_distribution<double> zDist(zMin, zMax);
      for (int i = 0; i < count; i++) {
        double r = radius * std::sqrt(uDist(rng));
        double theta = angleDist(rng);
        positions[i] =
            TVector3(r * std::cos(theta), r * std::sin(theta), zDist(rng));
      }
    } else {
      throw std::runtime_error(
          "electrons.position: expected [x,y,z], disk, or cylinder");
    }
  } else {
    throw std::runtime_error(
        "electrons.position: expected [x,y,z], disk, or cylinder");
  }

  return positions;
}

std::vector<ElectronConfig> GenerateElectrons(const YAML::Node& node) {
  if (!node["count"]) {
    throw std::runtime_error("electrons: missing required key 'count'");
  }
  int count = node["count"].as<int>();
  if (count <= 0) {
    throw std::runtime_error("electrons: count must be positive");
  }

  // Set up RNG
  std::mt19937 rng;
  if (node["seed"]) {
    rng.seed(node["seed"].as<unsigned int>());
  } else {
    std::random_device rd;
    rng.seed(rd());
  }

  // Resolve kinetic energy
  if (!node["kinetic_energy"]) {
    throw std::runtime_error(
        "electrons: missing required key 'kinetic_energy'");
  }
  if (!node["position"]) {
    throw std::runtime_error("electrons: missing required key 'position'");
  }

  bool hasPitch = node["pitch_angle"].IsDefined();
  bool hasIsotropic =
      node["isotropic"].IsDefined() && node["isotropic"].as<bool>();

  if (hasPitch && hasIsotropic) {
    throw std::runtime_error(
        "electrons: specify either 'pitch_angle' or 'isotropic: true', "
        "not both");
  }
  if (!hasPitch && !hasIsotropic) {
    throw std::runtime_error(
        "electrons: must specify either 'pitch_angle' or 'isotropic: true'");
  }

  auto keValues =
      ResolveScalarParam(node["kinetic_energy"], "kinetic_energy", count, rng);
  auto positions = ResolvePositions(node["position"], count, rng);

  // Build electron configs
  std::vector<ElectronConfig> electrons(count);
  if (hasPitch) {
    auto pitchValues =
        ResolveScalarParam(node["pitch_angle"], "pitch_angle", count, rng);
    for (int i = 0; i < count; i++) {
      electrons[i].position = positions[i];
      electrons[i].velocity =
          VelocityFromKEAndPitch(keValues[i], pitchValues[i]);
    }
  } else {
    // Isotropic: uniform random directions on the unit sphere
    std::uniform_real_distribution<double> cosThetaDist(-1.0, 1.0);
    std::uniform_real_distribution<double> phiDist(0.0, 2.0 * PI);
    for (int i = 0; i < count; i++) {
      double speed = GetSpeedFromKE(keValues[i], ME);
      double cosTheta = cosThetaDist(rng);
      double sinTheta = std::sqrt(1.0 - cosTheta * cosTheta);
      double phi = phiDist(rng);
      electrons[i].position = positions[i];
      electrons[i].velocity =
          TVector3(speed * sinTheta * std::cos(phi),
                   speed * sinTheta * std::sin(phi), speed * cosTheta);
    }
  }

  return electrons;
}

}  // namespace config
}  // namespace rad
