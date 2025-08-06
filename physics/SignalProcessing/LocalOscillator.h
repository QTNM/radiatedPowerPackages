/*
  LocalOscillator.h
  Simple model of a local oscillator used for downmixing of signals
*/

#ifndef LOCAL_OSCILLATOR_H
#define LOCAL_OSCILLATOR_H

#include "TMath.h"

namespace rad
{
  class LocalOscillator
  {
  private:
    double angularFreq;

  public:
    LocalOscillator();
    LocalOscillator(double angFreq);
    ~LocalOscillator();
    double GetAngularFrequency() { return angularFreq; }
    double GetFrequency() { return (angularFreq / (2*TMath::Pi())); }
    void SetAngularFrequency(double newFreq) { angularFreq = newFreq; }

    double GetInPhaseComponent(const double time);
    double GetQuadratureComponent(const double time);
  };
}

#endif
