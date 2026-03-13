#include "utilities/ConfigParser/ScatteringConfig.h"

#include <stdexcept>
#include <unordered_set>

#include "utilities/ConfigParser/YamlUtils.h"

namespace rad {
namespace config {

namespace {

const std::unordered_set<std::string> kValidSpecies = {
    "H", "H2", "He", "T", "T2"};

}  // namespace

ScatteringConfig ParseScatteringConfig(const YAML::Node& node) {
  ScatteringConfig config;

  if (!node["gases"] || !node["gases"].IsSequence()) {
    throw std::runtime_error(
        "scattering: missing required 'gases' sequence");
  }

  const YAML::Node& gases = node["gases"];
  if (gases.size() == 0) {
    throw std::runtime_error("scattering: 'gases' must contain at least one entry");
  }

  for (size_t i = 0; i < gases.size(); i++) {
    const auto& g = gases[i];
    std::string ctx = "scattering.gases[" + std::to_string(i) + "]";

    GasSpecies gs;
    gs.species = GetRequired<std::string>(g, "species", ctx);
    gs.density = GetRequired<double>(g, "density", ctx);

    if (kValidSpecies.find(gs.species) == kValidSpecies.end()) {
      throw std::runtime_error(ctx + ": unknown species '" + gs.species +
                               "' (valid: H, H2, He, T, T2)");
    }
    if (gs.density <= 0.0) {
      throw std::runtime_error(ctx + ": density must be positive");
    }

    config.gases.push_back(std::move(gs));
  }

  return config;
}

}  // namespace config
}  // namespace rad
