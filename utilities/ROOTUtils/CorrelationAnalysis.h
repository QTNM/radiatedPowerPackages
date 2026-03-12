/// CorrelationAnalysis.h - ROOT-based correlation analysis functions
#ifndef CORRELATION_ANALYSIS_H
#define CORRELATION_ANALYSIS_H

#include "TGraph.h"

namespace rad {

/// @brief Compute correlation of two arrays
/// @param length Length of both arrays
/// @param oldY1 First array
/// @param oldY2 Second array  
/// @return Correlation array of length *length*
double* GetCorrelation(int length, double* oldY1, double* oldY2);

/// @brief Compute correlation of two TGraphs
/// @param gr1 First input graph
/// @param gr2 Second input graph
/// @param zeroOffset Pointer to store zero-offset index (optional)
/// @return Correlation graph
TGraph* GetCorrelationGraph(const TGraph* gr1, const TGraph* gr2, int* zeroOffset = 0);

/// @brief Compute normalized correlation of two TGraphs
/// @param gr1 First input graph
/// @param gr2 Second input graph
/// @param zeroOffset Pointer to store zero-offset index (optional)
/// @return Normalized correlation graph
TGraph* GetNormalisedCorrelationGraph(const TGraph* gr1, const TGraph* gr2, int* zeroOffset = 0);

/// @brief Compute normalized correlation with time domain constraints
/// @param gr1 First input graph (must be zero-mean)
/// @param gr2 Second input graph (must be zero-mean)
/// @param zeroOffset Pointer to store zero-offset index (optional)
/// @param useDtRange Flag to enable delta-t range limits
/// @param dtMin Minimum delta-t to include
/// @param dtMax Maximum delta-t to include
/// @return Normalized correlation graph
TGraph* GetNormalisedCorrelationGraphTimeDomain(
    const TGraph* gr1, const TGraph* gr2, int* zeroOffset = 0,
    int useDtRange = 0, double dtMin = -1000, double dtMax = 1000);

}  // namespace rad

#endif