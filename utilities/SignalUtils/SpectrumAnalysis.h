/// SpectrumAnalysis.h - ROOT-free FFT and frequency domain analysis
#ifndef SPECTRUM_ANALYSIS_H
#define SPECTRUM_ANALYSIS_H

#include <vector>
#include <utility>

namespace rad {

/// Data structure for frequency domain results
/// First vector is frequencies [Hz], second is amplitudes/power
using FrequencySpectrum = std::pair<std::vector<double>, std::vector<double>>;

/// @brief Create power spectrum with time-integral normalization
/// @param timeData Time series data
/// @param deltaT Time step in seconds
/// @return Frequency spectrum pair (frequency [Hz], power)
FrequencySpectrum MakePowerSpectrumNorm(const std::vector<double>& timeData, double deltaT);

/// @brief Create power spectrum as periodogram
/// @param timeData Time series data
/// @param deltaT Time step in seconds
/// @return Frequency spectrum pair (frequency [Hz], power)
FrequencySpectrum MakePowerSpectrumPeriodogram(const std::vector<double>& timeData, double deltaT);

/// @brief Integrate power spectrum with normalization
/// @param freqData Frequency values [Hz]
/// @param powerData Power values
/// @param firstBin First bin to include (-1 for start)
/// @param lastBin Last bin to include (-1 for end)
/// @return Integrated power
double IntegratePowerNorm(const std::vector<double>& freqData, 
                         const std::vector<double>& powerData, 
                         int firstBin = -1, int lastBin = -1);

/// @brief Band pass filter using FFT
/// @param timeData Time series data
/// @param deltaT Time step in seconds
/// @param minFreq Lower cutoff frequency in Hz
/// @param maxFreq Upper cutoff frequency in Hz
/// @return Filtered time series data
std::vector<double> BandPassFilterSpectrum(const std::vector<double>& timeData, double deltaT,
                                          double minFreq, double maxFreq);

/// @brief Create FFT magnitude spectrum
/// @param timeData Time series data
/// @param deltaT Time step in seconds
/// @return Frequency spectrum pair (frequency [Hz], magnitude)
FrequencySpectrum MakeFFTMagnitude(const std::vector<double>& timeData, double deltaT);

/// @brief Add white noise to frequency domain power spectrum
/// @param powerData Power spectrum data (modified in place)
/// @param deltaF Frequency resolution in Hz
/// @param Teff Effective noise temperature in Kelvin
/// @param seed Random seed (0 for automatic)
void AddWhiteNoiseFrequencyDomain(std::vector<double>& powerData, double deltaF, 
                                 double Teff, int seed = 0);

}  // namespace rad

#endif  // SPECTRUM_ANALYSIS_H