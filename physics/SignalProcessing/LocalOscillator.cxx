// LocalOscillator.cxx

#include "physics/SignalProcessing/LocalOscillator.h"

#include "TMath.h"

rad::LocalOscillator::LocalOscillator(double angFreq) {
  angularFreq = angFreq;
}

rad::LocalOscillator::LocalOscillator() {
  angularFreq = 0;
}

rad::LocalOscillator::~LocalOscillator() { }

double rad::LocalOscillator::GetInPhaseComponent(const double time) {
  return (TMath::Cos(angularFreq * time));
}

double rad::LocalOscillator::GetQuadratureComponent(const double time) {
  return (TMath::Sin(angularFreq * time));  
}
