#ifndef DETECTOR_FACTORY_H
#define DETECTOR_FACTORY_H

#include <memory>
#include <variant>
#include <vector>
#include <yaml-cpp/yaml.h>

#include "physics/Antennas/IAntenna.h"
#include "physics/Waveguides/ICavity.h"
#include "physics/Waveguides/IWaveguide.h"
#include "physics/Waveguides/Probe.h"

namespace rad {
namespace config {

/// Antenna-based detector (one or more antennas)
struct AntennaDetector {
  std::vector<std::unique_ptr<IAntenna>> antennas;
};

/// Waveguide-based detector with readout probes
struct WaveguideDetector {
  std::unique_ptr<IWaveguide> waveguide;
  std::vector<Probe> probes;
};

/// Cavity-based detector
struct CavityDetector {
  std::unique_ptr<ICavity> cavity;
};

/// Variant unifying all detector types
using Detector =
    std::variant<AntennaDetector, WaveguideDetector, CavityDetector>;

/// @brief Create a Detector from a YAML detector configuration node
///
/// The node must contain a "type" key. Recognised types:
///   Antennas:   HalfWaveDipole, HertzianDipole, PatchAntenna, IsotropicAntenna
///   Waveguides: CircularWaveguide, RectangularWaveguide
///   Cavities:   CircularCavity
///
/// For waveguide types, a "probes" sequence is required.
///
/// @param detectorNode The YAML node containing "type" and "parameters"
/// @return Detector variant wrapping the constructed detector
/// @throws std::runtime_error if type is unknown or parameters are invalid
Detector CreateDetector(const YAML::Node& detectorNode);

}  // namespace config
}  // namespace rad

#endif
