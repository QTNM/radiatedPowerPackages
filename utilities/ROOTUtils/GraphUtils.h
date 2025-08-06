/// GraphUtils.h - ROOT TGraph utility functions
#ifndef GRAPH_UTILS_H
#define GRAPH_UTILS_H

#include <memory>
#include <vector>

#include "TGraph.h"
#include "TH1.h"

namespace rad {

/// @brief Set standard graph attributes for TGraph pointer
/// @param gr Graph to be formatted
void setGraphAttr(TGraph* gr);

/// @brief Set standard graph attributes for TGraph unique_ptr
/// @param gr Graph to be formatted  
void setGraphAttr(std::unique_ptr<TGraph>& gr);

/// @brief Set standard graph attributes for TGraph reference
/// @param gr Graph to be formatted
void setGraphAttr(TGraph& gr);

/// @brief Downmixes a time series with in-phase component
/// @param grInput Input time series graph
/// @param freq Downmix frequency in Hz
/// @return Downmixed graph
TGraph* DownmixInPhase(TGraph* grInput, double freq);

/// @brief Downmixes a time series with quadrature component  
/// @param grInput Input time series graph
/// @param freq Downmix frequency in Hz
/// @return Downmixed graph
TGraph* DownmixQuadrature(TGraph* grInput, double freq);

/// @brief Scales all points of a TGraph by a constant
/// @param grInput The input graph
/// @param scale Scale factor
void ScaleGraph(TGraph* grInput, double scale);

/// @brief Sums multiple TGraphs point by point
/// @param grInput Vector of graphs to be summed
/// @return Summed graph
TGraph* SumGraphs(std::vector<TGraph*> grInput);

/// @brief Samples waveform at given rate using linear interpolation
/// @param grInput Input graph to be sampled
/// @param sRate Sample rate in Hz
/// @return Downsampled graph
TGraph* SampleWaveform(TGraph* grInput, double sRate);

/// @brief Converts TGraph to histogram
/// @param grInput Input graph to be converted
/// @return Converted histogram
TH1D* GraphToHistogram(TGraph* grInput);

/// @brief Processes signal: downmix, filter, and sample
/// @param grInput Input time domain voltage graph
/// @param downmixFreq Downmix frequency in Hz
/// @param sampleRate Sample rate in Hz
/// @return Processed time domain graph
TGraph* SignalProcessGraph(TGraph* grInput, double downmixFreq, double sampleRate);

/// @brief Sum power in a power spectrum graph
/// @param gr Pointer to input power spectrum TGraph
/// @param firstBin First bin to include (-1 for start)
/// @param lastBin Last bin to include (-1 for end)
/// @return Sum of power values
double SumPower(const TGraph* gr, int firstBin = -1, int lastBin = -1);

/// @brief Sum voltage squared in a time series
/// @param gr Pointer to input TGraph
/// @param firstBin First bin to include (-1 for start)
/// @param lastBin Last bin to include (-1 for end)
/// @return Sum of voltage squared
double SumVoltageSquared(const TGraph* gr, int firstBin = -1, int lastBin = -1);

}  // namespace rad

#endif