/*
  ScaledSignal.h

  Derived class of Signal which allows for abitrary scaling of the input voltage signal
*/

#ifndef SCALED_SIGNAL_H
#define SCALED_SIGNAL_H

#include "physics/SignalProcessing/Signal.h"
#include "physics/SignalProcessing/LocalOscillator.h"
#include "physics/SignalProcessing/NoiseFunc.h"
#include "physics/SignalProcessing/InducedVoltage.h"

#include <vector>

namespace rad
{
  class ScaledSignal : public Signal {

  private:
    double scaleFactor;

    void ProcessTimeChunk(InducedVoltage iv, LocalOscillator lo, double thisChunk, double lastChunk,
			  std::vector<GaussianNoise> noiseTerms,
			  double &firstSampleTime, double &firstSample10Time);
        
  public:
    /// Constructor for a single induced voltage
    /// \param iv The induced voltage at a particular field point
    /// \param lo The local oscillator used to do the downmixing
    /// \param srate The sample rate of a hypothetical ADC
    /// \param scaleFactor Factor by which to multiply the original signal
    /// \param noiseTerms A vector of noise terms to be added to the signal
    ScaledSignal(InducedVoltage iv, LocalOscillator lo, double srate, double scaleFactor,
		 std::vector<GaussianNoise> noiseTerms={}, double maxTime=-1);

    /// Copy constructor
    ScaledSignal(const ScaledSignal &s1);
  };
}

#endif
