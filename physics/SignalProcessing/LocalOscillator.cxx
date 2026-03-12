// LocalOscillator.cxx

#include "physics/SignalProcessing/LocalOscillator.h"

#include <cmath>

rad::LocalOscillator::LocalOscillator(double angFreq) {
  angularFreq = angFreq;
}

rad::LocalOscillator::LocalOscillator() {
  angularFreq = 0;
}

rad::LocalOscillator::~LocalOscillator() { }

double rad::LocalOscillator::GetInPhaseComponent(const double time) {
  return (std::cos(angularFreq * time));
}

double rad::LocalOscillator::GetQuadratureComponent(const double time) {
  return (std::sin(angularFreq * time));  
}
