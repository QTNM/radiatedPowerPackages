#include "Waveguides/Probe.h"

rad::Probe::Probe(TVector3 pos, WaveguideMode m, bool state)
    : position(pos), mode(m), polarisationState(state) {}

rad::Probe::Probe()
    : position(TVector3(0, 0, 0)),
      mode(WaveguideMode()),
      polarisationState(true) {}

TVector3 rad::Probe::GetPosition() { return position; }

rad::WaveguideMode rad::Probe::GetMode() { return mode; }

bool rad::Probe::GetPolarisationState() { return polarisationState; }