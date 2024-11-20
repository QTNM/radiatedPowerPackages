#include "Waveguides/Probe.h"

rad::Probe::Probe(TVector3 pos, WaveguideMode m, bool state)
    : position(pos), mode(m), polarisationState(state) {}

rad::Probe::Probe()
    : position(TVector3(0, 0, 0)),
      mode(WaveguideMode()),
      polarisationState(true) {}

inline TVector3 rad::Probe::GetPosition() { return position; }

inline rad::WaveguideMode rad::Probe::GetMode() { return mode; }

inline bool rad::Probe::GetPolarisationState() { return polarisationState; }