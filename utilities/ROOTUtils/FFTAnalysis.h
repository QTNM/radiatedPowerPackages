/// FFTAnalysis.h - ROOT-based FFT and frequency domain analysis
#ifndef FFT_ANALYSIS_H
#define FFT_ANALYSIS_H

#include "TGraph.h"

namespace rad {

/// @brief Create power spectrum with time-integral normalization  
/// @param grWave Input time series graph
/// @return Power spectrum graph
TGraph* MakePowerSpectrumNorm(const TGraph* grWave);

/// @brief Create power spectrum as periodogram
/// @param grWave Input time series graph
/// @return Periodogram graph (pointer)
TGraph* MakePowerSpectrumPeriodogram(const TGraph* grWave);

/// @brief Create power spectrum as periodogram (by reference)
/// @param grWave Input time series graph
/// @return Periodogram graph (by value)
TGraph MakePowerSpectrumPeriodogram(const TGraph& grWave);

/// @brief Integrate power spectrum with normalization
/// @param grFFT Power spectrum graph
/// @param firstBin First bin to include (-1 for start)
/// @param lastBin Last bin to include (-1 for end)
/// @return Integrated power
double IntegratePowerNorm(const TGraph* grFFT, Int_t firstBin = -1, Int_t lastBin = -1);

/// @brief Band pass filter using FFT
/// @param grWave Input time series graph
/// @param minFreq Lower cutoff frequency in Hz
/// @param maxFreq Upper cutoff frequency in Hz
/// @return Filtered graph
TGraph* BandPassFilter(const TGraph* grWave, double minFreq, double maxFreq);

/// @brief Create FFT magnitude spectrum
/// @param grInput Input time series graph
/// @return FFT magnitude graph
TGraph* MakeFFTMagGraph(TGraph* grInput);

/// @brief Add white noise to frequency domain power spectrum
/// @param grIn Input power spectrum graph (modified in place)
/// @param Teff Effective noise temperature in Kelvin
/// @param seed Random seed (0 for automatic)
void AddWhiteNoiseFrequencyDomainPowerNorm(TGraph* grIn, double Teff, int seed = 0);

}  // namespace rad

#endif