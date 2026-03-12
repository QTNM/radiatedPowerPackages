/// SpectrumAnalysis.cxx - ROOT-free FFT and frequency domain analysis
#include "utilities/SignalUtils/SpectrumAnalysis.h"

#include "utilities/BasicCore/Constants.h"
#include "utilities/SignalUtils/FourierTransforms.h"
#include "utilities/SignalUtils/FFTWComplex.h"
#include <cmath>
#include <random>
#include <algorithm>

using namespace rad;

FrequencySpectrum rad::MakePowerSpectrumNorm(const std::vector<double>& timeData, double deltaT) {
    const int length = timeData.size();
    const double* dataPtr = timeData.data();
    
    FFTWComplex* theFFT = doFFT(length, const_cast<double*>(dataPtr));
    
    const int newLength = (length / 2) + 1;
    std::vector<double> frequencies(newLength);
    std::vector<double> powers(newLength);
    
    const double deltaF = 1.0 / (deltaT * length);
    
    for (int i = 0; i < newLength; i++) {
        double power = std::pow(getAbs(theFFT[i]), 2);
        if (i > 0 && i < newLength - 1) {
            power *= 2;  // account for symmetry
        }
        power *= deltaT / length;  // For time-integral squared amplitude
        power /= deltaF;           // Just to normalise bin-widths
        
        frequencies[i] = i * deltaF;
        powers[i] = power;
    }
    
    delete[] theFFT;
    return std::make_pair(frequencies, powers);
}

FrequencySpectrum rad::MakePowerSpectrumPeriodogram(const std::vector<double>& timeData, double deltaT) {
    const int length = timeData.size();
    const double* dataPtr = timeData.data();
    
    FFTWComplex* theFFT = doFFT(length, const_cast<double*>(dataPtr));
    
    const int newLength = (length / 2) + 1;
    std::vector<double> frequencies(newLength);
    std::vector<double> powers(newLength);
    
    const double deltaF = 1.0 / (deltaT * length);
    const double lengthDub = static_cast<double>(length);
    
    for (int i = 0; i < newLength; i++) {
        double power = std::pow(getAbs(theFFT[i]), 2);
        if (i > 0 && i < newLength - 1) {
            power *= 2;  // account for symmetry
        }
        const double scale = lengthDub * lengthDub;
        power /= scale;
        
        frequencies[i] = i * deltaF;
        powers[i] = power;
    }
    
    delete[] theFFT;
    return std::make_pair(frequencies, powers);
}

double rad::IntegratePowerNorm(const std::vector<double>& freqData, 
                              const std::vector<double>& powerData, 
                              int firstBin, int lastBin) {
    const int nPoints = static_cast<int>(powerData.size());
    
    // Handle default values
    if (firstBin < 0) firstBin = 0;
    if (lastBin < 0) lastBin = nPoints - 1;
    
    // Ensure valid range
    firstBin = std::max(0, firstBin);
    lastBin = std::min(nPoints - 1, lastBin);
    
    if (firstBin >= lastBin) return 0.0;
    
    double integral = 0.0;
    const double deltaF = (freqData.size() > 1) ? (freqData[1] - freqData[0]) : 1.0;
    
    for (int i = firstBin; i <= lastBin; i++) {
        integral += powerData[i] * deltaF;
    }
    
    integral *= deltaF;  // Additional factor as in original
    return integral;
}

std::vector<double> rad::BandPassFilterSpectrum(const std::vector<double>& timeData, double deltaT,
                                               double minFreq, double maxFreq) {
    const int length = timeData.size();
    const double* dataPtr = timeData.data();
    
    FFTWComplex* theFFT = doFFT(length, const_cast<double*>(dataPtr));
    
    const int newLength = (length / 2) + 1;
    const double deltaF = 1.0 / (deltaT * length);
    
    // Apply filter in frequency domain
    for (int i = 0; i < newLength; i++) {
        const double freq = i * deltaF;
        if (freq < minFreq || freq > maxFreq) {
            theFFT[i].re = 0.0;
            theFFT[i].im = 0.0;
        }
    }
    
    double* filteredData = doInverseFFT(length, theFFT);
    
    std::vector<double> result(filteredData, filteredData + length);
    
    delete[] theFFT;
    delete[] filteredData;
    
    return result;
}

FrequencySpectrum rad::MakeFFTMagnitude(const std::vector<double>& timeData, double deltaT) {
    const int length = timeData.size();
    const double* dataPtr = timeData.data();
    
    FFTWComplex* theFFT = doFFT(length, const_cast<double*>(dataPtr));
    
    const int newLength = (length / 2) + 1;
    std::vector<double> frequencies(newLength);
    std::vector<double> magnitudes(newLength);
    
    const double deltaF = 1.0 / (deltaT * length);
    
    for (int i = 0; i < newLength; i++) {
        frequencies[i] = i * deltaF;
        magnitudes[i] = getAbs(theFFT[i]);
    }
    
    delete[] theFFT;
    return std::make_pair(frequencies, magnitudes);
}

void rad::AddWhiteNoiseFrequencyDomain(std::vector<double>& powerData, double deltaF, 
                                      double Teff, int seed) {
    // Set up random number generation
    std::mt19937 gen;
    if (seed != 0) {
        gen.seed(seed);
    } else {
        std::random_device rd;
        gen.seed(rd());
    }
    
    // Calculate noise parameters
    const double sampleRate = 2.0 * (powerData.size() - 1) * deltaF;
    const double deltaT = 1.0 / sampleRate;
    const double sigma = std::sqrt(K_B * Teff * sampleRate);
    
    // Rayleigh distribution
    std::exponential_distribution<double> exp_dist(1.0 / (sigma * sigma));
    
    // Add noise to each bin
    for (size_t i = 0; i < powerData.size(); i++) {
        // Generate Rayleigh distributed noise using exponential distribution
        const double u = exp_dist(gen);
        const double rayleigh_sample = std::sqrt(2.0 * u) * sigma;
        
        const double noise = std::pow(rayleigh_sample * std::sqrt(0.5), 2) * (1.0 / deltaF) * deltaT;
        powerData[i] += noise;
    }
}