#ifndef FIELD_FACTORY_H
#define FIELD_FACTORY_H

#include <memory>
#include <yaml-cpp/yaml.h>

#include "physics/ElectronDynamics/BaseField.h"

namespace rad {
namespace config {

/// @brief Create a BaseField from a YAML field configuration node
///
/// The node must contain a "type" key naming the field class, and a
/// "parameters" node whose keys match the C++ constructor parameter names.
///
/// Currently supported types: UniformField, BathtubField, HarmonicField
///
/// @param fieldNode The YAML node containing "type" and "parameters"
/// @return Unique pointer to the constructed field
/// @throws std::runtime_error if type is unknown or parameters are invalid
std::unique_ptr<BaseField> CreateField(const YAML::Node& fieldNode);

}  // namespace config
}  // namespace rad

#endif
