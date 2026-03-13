#include "utilities/ConfigParser/DetectorFactory.h"

#include <functional>
#include <memory>
#include <stdexcept>
#include <string>
#include <unordered_map>

#include "physics/Antennas/HalfWaveDipole.h"
#include "physics/Antennas/HertzianDipole.h"
#include "physics/Antennas/IsotropicAntenna.h"
#include "physics/Antennas/PatchAntenna.h"
#include "physics/Waveguides/CircularCavity.h"
#include "physics/Waveguides/CircularWaveguide.h"
#include "physics/Waveguides/RectangularWaveguide.h"
#include "physics/Waveguides/WaveguideMode.h"
#include "utilities/ConfigParser/YamlUtils.h"

namespace rad {
namespace config {

// --- Antenna registry ---

using AntennaCreator =
    std::function<std::unique_ptr<IAntenna>(const YAML::Node&)>;

static const std::unordered_map<std::string, AntennaCreator>&
GetAntennaRegistry() {
  static const std::unordered_map<std::string, AntennaCreator> registry = {
      {"HalfWaveDipole",
       [](const YAML::Node& params) {
         TVector3 pos =
             ParseVector3(params["position"], "HalfWaveDipole.position");
         TVector3 xAxis =
             ParseVector3(params["x_axis"], "HalfWaveDipole.x_axis");
         TVector3 zAxis =
             ParseVector3(params["z_axis"], "HalfWaveDipole.z_axis");
         double freq =
             GetRequired<double>(params, "frequency", "HalfWaveDipole");
         double delay = GetOptional<double>(params, "delay", 0.0);
         return std::make_unique<HalfWaveDipole>(pos, xAxis, zAxis, freq,
                                                  delay);
       }},

      {"HertzianDipole",
       [](const YAML::Node& params) {
         TVector3 pos =
             ParseVector3(params["position"], "HertzianDipole.position");
         TVector3 xAxis =
             ParseVector3(params["x_axis"], "HertzianDipole.x_axis");
         TVector3 zAxis =
             ParseVector3(params["z_axis"], "HertzianDipole.z_axis");
         double freq =
             GetRequired<double>(params, "frequency", "HertzianDipole");
         double delay = GetOptional<double>(params, "delay", 0.0);
         return std::make_unique<HertzianDipole>(pos, xAxis, zAxis, freq,
                                                  delay);
       }},

      {"PatchAntenna",
       [](const YAML::Node& params) {
         TVector3 pos =
             ParseVector3(params["position"], "PatchAntenna.position");
         TVector3 xAxis =
             ParseVector3(params["x_axis"], "PatchAntenna.x_axis");
         TVector3 yAxis =
             ParseVector3(params["y_axis"], "PatchAntenna.y_axis");
         double width =
             GetRequired<double>(params, "width", "PatchAntenna");
         double height =
             GetRequired<double>(params, "height", "PatchAntenna");
         double freq =
             GetRequired<double>(params, "frequency", "PatchAntenna");
         double perm = GetOptional<double>(params, "permittivity", 1.0);
         double delay = GetOptional<double>(params, "delay", 0.0);
         return std::make_unique<PatchAntenna>(pos, xAxis, yAxis, width,
                                                height, freq, perm, delay);
       }},

      {"IsotropicAntenna",
       [](const YAML::Node& params) {
         TVector3 pos =
             ParseVector3(params["position"], "IsotropicAntenna.position");
         double aEff =
             GetRequired<double>(params, "effective_area", "IsotropicAntenna");
         double Z =
             GetRequired<double>(params, "impedance", "IsotropicAntenna");
         double freq = GetOptional<double>(params, "frequency", 27e9);
         double delay = GetOptional<double>(params, "delay", 0.0);
         return std::make_unique<IsotropicAntenna>(pos, aEff, Z, freq, delay);
       }},
  };
  return registry;
}

// --- Waveguide registry ---

using WaveguideCreator =
    std::function<std::unique_ptr<IWaveguide>(const YAML::Node&)>;

static const std::unordered_map<std::string, WaveguideCreator>&
GetWaveguideRegistry() {
  static const std::unordered_map<std::string, WaveguideCreator> registry = {
      {"CircularWaveguide",
       [](const YAML::Node& params) {
         double radius =
             GetRequired<double>(params, "radius", "CircularWaveguide");
         double length =
             GetRequired<double>(params, "length", "CircularWaveguide");
         return std::make_unique<CircularWaveguide>(radius, length);
       }},

      {"RectangularWaveguide",
       [](const YAML::Node& params) {
         double longSide =
             GetRequired<double>(params, "long_side", "RectangularWaveguide");
         double shortSide =
             GetRequired<double>(params, "short_side", "RectangularWaveguide");
         double length =
             GetRequired<double>(params, "length", "RectangularWaveguide");
         return std::make_unique<RectangularWaveguide>(longSide, shortSide,
                                                        length);
       }},
  };
  return registry;
}

// --- Cavity registry ---

using CavityCreator =
    std::function<std::unique_ptr<ICavity>(const YAML::Node&)>;

static const std::unordered_map<std::string, CavityCreator>&
GetCavityRegistry() {
  static const std::unordered_map<std::string, CavityCreator> registry = {
      {"CircularCavity",
       [](const YAML::Node& params) {
         double radius =
             GetRequired<double>(params, "radius", "CircularCavity");
         double length =
             GetRequired<double>(params, "length", "CircularCavity");
         TVector3 probe = ParseVector3(params["probe_position"],
                                        "CircularCavity.probe_position");
         return std::make_unique<CircularCavity>(radius, length, probe);
       }},
  };
  return registry;
}

// --- Probe parsing helpers ---

static WaveguideMode ParseWaveguideMode(const YAML::Node& node,
                                         const std::string& context) {
  std::string typeStr = GetRequired<std::string>(node, "type", context);
  ModeType mt;
  if (typeStr == "TE")
    mt = kTE;
  else if (typeStr == "TM")
    mt = kTM;
  else if (typeStr == "TEM")
    mt = kTEM;
  else
    throw std::runtime_error(context + ": unknown mode type '" + typeStr + "'");

  unsigned int m = GetRequired<unsigned int>(node, "m", context);
  unsigned int n = GetRequired<unsigned int>(node, "n", context);
  return WaveguideMode(m, n, mt);
}

static Probe ParseProbe(const YAML::Node& node, const std::string& context) {
  TVector3 pos = ParseVector3(node["position"], context + ".position");
  if (!node["mode"]) {
    throw std::runtime_error(context + ": missing required key 'mode'");
  }
  WaveguideMode mode = ParseWaveguideMode(node["mode"], context + ".mode");
  bool state = GetOptional<bool>(node, "polarisation", true);
  return Probe(pos, mode, state);
}

// --- Public factory function ---

Detector CreateDetector(const YAML::Node& detectorNode) {
  if (!detectorNode["type"]) {
    throw std::runtime_error("detector: missing required key 'type'");
  }
  const std::string typeName = detectorNode["type"].as<std::string>();

  const YAML::Node& params = detectorNode["parameters"];
  if (!params) {
    throw std::runtime_error("detector: type '" + typeName +
                             "' requires a 'parameters' section");
  }

  // Check antenna registry first
  {
    const auto& registry = GetAntennaRegistry();
    auto it = registry.find(typeName);
    if (it != registry.end()) {
      AntennaDetector det;
      det.antennas.push_back(it->second(params));
      return det;
    }
  }

  // Check waveguide registry
  {
    const auto& registry = GetWaveguideRegistry();
    auto it = registry.find(typeName);
    if (it != registry.end()) {
      WaveguideDetector det;
      det.waveguide = it->second(params);

      if (!detectorNode["probes"] || !detectorNode["probes"].IsSequence() ||
          detectorNode["probes"].size() == 0) {
        throw std::runtime_error(
            "detector: waveguide type '" + typeName +
            "' requires a non-empty 'probes' sequence");
      }
      for (size_t i = 0; i < detectorNode["probes"].size(); i++) {
        det.probes.push_back(
            ParseProbe(detectorNode["probes"][i],
                        "detector.probes[" + std::to_string(i) + "]"));
      }
      return det;
    }
  }

  // Check cavity registry
  {
    const auto& registry = GetCavityRegistry();
    auto it = registry.find(typeName);
    if (it != registry.end()) {
      CavityDetector det;
      det.cavity = it->second(params);
      return det;
    }
  }

  throw std::runtime_error("detector: unknown type '" + typeName + "'");
}

}  // namespace config
}  // namespace rad
