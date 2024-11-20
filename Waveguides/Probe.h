/*
  Probe.h

  Class containing information about probes for use in waveguide simulations
*/

#ifndef PROBE_H
#define PROBE_H

#include "TVector3.h"
#include "Waveguides/WaveguideMode.h"

namespace rad {
class Probe {
 private:
  TVector3 position;       // Probe position
  WaveguideMode mode;      // Mode read out by probe
  bool polarisationState;  // Polarisation state of probe (where appropriate)

 public:
  /// @brief Parametrised constructor
  /// @param pos Three vector for probe position (metres)
  /// @param m Mode which the probe reads out
  /// @param state Polarisation state of the probe (where applicable
  Probe(TVector3 pos, WaveguideMode m, bool state = true);

  /// @brief Default constructor
  Probe();

  /// @brief Getter for probe position
  /// @return The position of the probe
  TVector3 GetPosition();

  /// @brief Getter for probe mode
  /// @return The mode of the probe
  WaveguideMode GetMode();

  /// @brief Getter for probe polarisation state
  /// @return The polarisation state of the probe
  bool GetPolarisationState();
};
}  // namespace rad

#endif