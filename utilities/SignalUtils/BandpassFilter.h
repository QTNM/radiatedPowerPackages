/// BandpassFilter.h
#ifndef BANDPASS_FILTER_H
#define BANDPASS_FILTER_H

#include <iostream>
#include <vector>

#include "utilities/SignalUtils/FFTWComplex.h"
#include "utilities/SignalUtils/FourierTransforms.h"

namespace rad {
/// @brief Implements band pass filter on vectors of values
/// @param xVals Vector of equally spaced time values
/// @param yVals Vector of corresponding y values
/// @param minFreq Lower cutoff frequency in Hz
/// @param maxFreq Upper cutoff frequency in Hz
/// @return Vector of filtered y values
std::vector<double> BandPassFilter(std::vector<double> xVals,
                                   std::vector<double> yVals, double minFreq,
                                   double maxFreq);
}  // namespace rad

#endif