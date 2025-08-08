/// FFTAnalysis.cxx - ROOT-based FFT and frequency domain analysis (using SignalUtils internally)
#include "utilities/ROOTUtils/FFTAnalysis.h"

#include "utilities/BasicCore/Constants.h"
#include "utilities/ROOTUtils/GraphUtils.h"
#include "utilities/SignalUtils/SpectrumAnalysis.h"
#include <cmath>
#include "TF1.h"
#include "TRandom3.h"

using namespace rad;

TGraph *rad::MakePowerSpectrumNorm(const TGraph *grWave) {
  // Extract data from TGraph
  double *timeData = grWave->GetY();
  double *xData = grWave->GetX();
  int n = grWave->GetN();
  double deltaT = xData[1] - xData[0];
  
  // Convert to vector and process using SignalUtils
  std::vector<double> timeVec(timeData, timeData + n);
  auto spectrum = rad::MakePowerSpectrumNorm(timeVec, deltaT);
  
  // Convert back to TGraph
  TGraph *result = new TGraph(spectrum.first.size(), spectrum.first.data(), spectrum.second.data());
  setGraphAttr(result);
  return result;
}

TGraph *rad::MakePowerSpectrumPeriodogram(const TGraph *grWave) {
  // Extract data from TGraph
  double *timeData = grWave->GetY();
  double *xData = grWave->GetX();
  int n = grWave->GetN();
  double deltaT = xData[1] - xData[0];
  
  // Convert to vector and process using SignalUtils
  std::vector<double> timeVec(timeData, timeData + n);
  auto spectrum = rad::MakePowerSpectrumPeriodogram(timeVec, deltaT);
  
  // Convert back to TGraph
  TGraph *result = new TGraph(spectrum.first.size(), spectrum.first.data(), spectrum.second.data());
  setGraphAttr(result);
  result->GetXaxis()->SetTitle("Frequency [Hz]");
  return result;
}

TGraph rad::MakePowerSpectrumPeriodogram(const TGraph &grWave) {
  // Extract data from TGraph
  double *timeData = grWave.GetY();
  double *xData = grWave.GetX();
  int n = grWave.GetN();
  double deltaT = xData[1] - xData[0];
  
  // Convert to vector and process using SignalUtils
  std::vector<double> timeVec(timeData, timeData + n);
  auto spectrum = rad::MakePowerSpectrumPeriodogram(timeVec, deltaT);
  
  // Convert back to TGraph
  TGraph result(spectrum.first.size(), spectrum.first.data(), spectrum.second.data());
  setGraphAttr(result);
  result.GetXaxis()->SetTitle("Frequency [Hz]");
  return result;
}

double rad::IntegratePowerNorm(const TGraph *grFFT, Int_t firstBin,
                               Int_t lastBin) {
  // Extract data from TGraph
  double *freqData = grFFT->GetX();
  double *powerData = grFFT->GetY();
  int n = grFFT->GetN();
  
  // Handle default values
  if (firstBin < 0) firstBin = 0;
  if (lastBin < 0) lastBin = n - 1;
  
  // Convert to vectors
  std::vector<double> freqVec(freqData, freqData + n);
  std::vector<double> powerVec(powerData, powerData + n);
  
  return rad::IntegratePowerNorm(freqVec, powerVec, firstBin, lastBin);
}

TGraph *rad::BandPassFilter(const TGraph *grWave, const double minFreq,
                            const double maxFreq) {
  // Extract data from TGraph
  double *timeData = grWave->GetY();
  double *xData = grWave->GetX();
  int n = grWave->GetN();
  double deltaT = xData[1] - xData[0];
  
  // Convert to vector and process using SignalUtils
  std::vector<double> timeVec(timeData, timeData + n);
  auto filteredData = rad::BandPassFilterSpectrum(timeVec, deltaT, minFreq, maxFreq);
  
  // Convert back to TGraph
  return new TGraph(n, xData, filteredData.data());
}

TGraph *rad::MakeFFTMagGraph(TGraph *grInput) {
  // Extract data from TGraph
  double *timeData = grInput->GetY();
  double *xData = grInput->GetX();
  int n = grInput->GetN();
  double deltaT = xData[1] - xData[0];
  
  // Convert to vector and process using SignalUtils
  std::vector<double> timeVec(timeData, timeData + n);
  auto spectrum = rad::MakeFFTMagnitude(timeVec, deltaT);
  
  // Convert back to TGraph
  TGraph *result = new TGraph(spectrum.first.size(), spectrum.first.data(), spectrum.second.data());
  setGraphAttr(result);
  return result;
}

/// PDF for Rayleigh distribution
double RayleighPDFFunc(double *x, double *par) {
  double xx = x[0];
  double sigma = par[0];
  double sigmaSq = sigma * sigma;
  double prob = (xx / sigmaSq) * exp(-xx * xx / (2 * sigmaSq));
  return prob;
}

void rad::AddWhiteNoiseFrequencyDomainPowerNorm(TGraph *grIn, const double Teff,
                                                const int seed) {
  // Extract data from TGraph
  double *powerData = grIn->GetY();
  int n = grIn->GetN();
  double deltaF = grIn->GetPointX(1) - grIn->GetPointX(0);
  
  // Convert to vector and process using SignalUtils
  std::vector<double> powerVec(powerData, powerData + n);
  rad::AddWhiteNoiseFrequencyDomain(powerVec, deltaF, Teff, seed);
  
  // Update TGraph with modified data
  for (int i = 0; i < n; i++) {
    grIn->SetPointY(i, powerVec[i]);
  }
}