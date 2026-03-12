#include "utilities/ConfigParser/FieldFactory.h"

#include <functional>
#include <memory>
#include <stdexcept>
#include <string>
#include <unordered_map>

#include "physics/ElectronDynamics/QTNMFields.h"
#include "utilities/ConfigParser/YamlUtils.h"

namespace rad {
namespace config {

using FieldCreator =
    std::function<std::unique_ptr<BaseField>(const YAML::Node&)>;

static const std::unordered_map<std::string, FieldCreator>& GetRegistry() {
  static const std::unordered_map<std::string, FieldCreator> registry = {
      {"UniformField",
       [](const YAML::Node& params) {
         double field = GetRequired<double>(params, "field", "UniformField");
         return std::make_unique<UniformField>(field);
       }},

      {"BathtubField",
       [](const YAML::Node& params) {
         double radius = GetRequired<double>(params, "radius", "BathtubField");
         double current =
             GetRequired<double>(params, "current", "BathtubField");
         double Z1 = GetRequired<double>(params, "Z1", "BathtubField");
         double Z2 = GetRequired<double>(params, "Z2", "BathtubField");
         TVector3 background =
             ParseVector3(params["background"], "BathtubField.background");
         return std::make_unique<BathtubField>(radius, current, Z1, Z2,
                                               background);
       }},

      {"HarmonicField",
       [](const YAML::Node& params) {
         double radius =
             GetRequired<double>(params, "radius", "HarmonicField");
         double current =
             GetRequired<double>(params, "current", "HarmonicField");
         double background =
             GetRequired<double>(params, "background", "HarmonicField");
         return std::make_unique<HarmonicField>(radius, current, background);
       }},
  };
  return registry;
}

std::unique_ptr<BaseField> CreateField(const YAML::Node& fieldNode) {
  if (!fieldNode["type"]) {
    throw std::runtime_error("field: missing required key 'type'");
  }
  const std::string typeName = fieldNode["type"].as<std::string>();

  const auto& registry = GetRegistry();
  auto it = registry.find(typeName);
  if (it == registry.end()) {
    throw std::runtime_error("field: unknown type '" + typeName + "'");
  }

  const YAML::Node& params = fieldNode["parameters"];
  if (!params) {
    throw std::runtime_error("field: type '" + typeName +
                             "' requires a 'parameters' section");
  }

  return it->second(params);
}

}  // namespace config
}  // namespace rad
