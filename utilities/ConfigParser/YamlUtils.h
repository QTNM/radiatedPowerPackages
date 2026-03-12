#ifndef YAML_UTILS_H
#define YAML_UTILS_H

#include <stdexcept>
#include <string>
#include <yaml-cpp/yaml.h>
#include "TVector3.h"

namespace rad {
namespace config {

/// @brief Parse a TVector3 from a YAML 3-element sequence [x, y, z]
/// @param node The YAML node (must be a sequence of 3 doubles)
/// @param context Description of where this vector is in the config, for errors
/// @return The parsed TVector3
/// @throws std::runtime_error if node is not a 3-element sequence
inline TVector3 ParseVector3(const YAML::Node& node,
                             const std::string& context) {
  if (!node || !node.IsSequence() || node.size() != 3) {
    throw std::runtime_error(context + ": expected a 3-element array [x, y, z]");
  }
  return TVector3(node[0].as<double>(), node[1].as<double>(),
                  node[2].as<double>());
}

/// @brief Get a required parameter from a YAML node
/// @tparam T The value type to extract
/// @param params The YAML node containing the parameter
/// @param key The parameter name
/// @param typeName The parent type name, for error messages
/// @return The parameter value
/// @throws std::runtime_error if the key is missing
template <typename T>
T GetRequired(const YAML::Node& params, const std::string& key,
              const std::string& typeName) {
  if (!params[key]) {
    throw std::runtime_error(typeName + ": missing required parameter '" + key +
                             "'");
  }
  return params[key].as<T>();
}

/// @brief Get an optional parameter from a YAML node, with a default value
/// @tparam T The value type to extract
/// @param params The YAML node containing the parameter
/// @param key The parameter name
/// @param defaultVal The value to return if the key is absent
/// @return The parameter value, or defaultVal if absent
template <typename T>
T GetOptional(const YAML::Node& params, const std::string& key, T defaultVal) {
  if (params[key]) {
    return params[key].as<T>();
  }
  return defaultVal;
}

}  // namespace config
}  // namespace rad

#endif
