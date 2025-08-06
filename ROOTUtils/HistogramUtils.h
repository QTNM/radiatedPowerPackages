/// HistogramUtils.h - ROOT histogram utility functions
#ifndef HISTOGRAM_UTILS_H
#define HISTOGRAM_UTILS_H

#include <memory>

#include "TH1.h"
#include "TH2.h"

namespace rad {

/// @brief Set standard histogram attributes for TH1 pointer
/// @param h Histogram to be formatted
void SetHistAttr(TH1* h);

/// @brief Set standard histogram attributes for TH1 reference
/// @param h Histogram to be formatted
void SetHistAttr(TH1& h);

/// @brief Set standard histogram attributes for TH1D unique_ptr
/// @param h Histogram to be formatted
void SetHistAttr(std::unique_ptr<TH1D>& h);

/// @brief Set standard histogram attributes for TH2 pointer
/// @param h Histogram to be formatted
void SetHistAttr(TH2* h);

/// @brief Set standard histogram attributes for TH2 unique_ptr
/// @param h Histogram to be formatted
void SetHistAttr(std::unique_ptr<TH2>& h);

}  // namespace rad

#endif