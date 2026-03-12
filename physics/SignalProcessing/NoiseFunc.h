/*
  NoiseFunc.h
  Contains functions for implementing noise in signal processing
*/
#ifndef NOISE_FUNC_H
#define NOISE_FUNC_H

#include "TRandom3.h"

namespace rad
{

  /// Gaussian random noise which is constant in time
  class GaussianNoise
  {
  private:
    TRandom3 *numGen = 0; 
    double noiseTemp;
    double sampleFreq;
    double resistance;
    double sigma;
    
    void SetResistance(double r);
    
  public:
    GaussianNoise();
    /// \param T is the noise temperature in kelvin
    /// \param R is resistance of the load circuit
    /// \param setSeed is the seed for the random number generator
    GaussianNoise(double T, double R, int setSeed=1234); 
    ~GaussianNoise();
    GaussianNoise(const GaussianNoise &g1);

    void SetSampleFreq(double fs);
    void SetSigma();
    double GetNoiseVoltage(bool IsComponent=true);
    double GetSigma() { return sigma; }
    double GetNoiseTemp() { return noiseTemp; }
    double GetFs() { return sampleFreq; }
    double GetResistance() { return resistance; }
  };
}

#endif
