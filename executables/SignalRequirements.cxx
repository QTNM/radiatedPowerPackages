/*
  SignalRequirements.cxx

  Using known noise properties, calculate the required signal power
  Allow for some hypothesised noise rate
*/

#include <iostream>
#include <memory>

#include "BasicFunctions/BasicFunctions.h"
#include "BasicFunctions/Constants.h"
#include "TF1.h"
#include "TFile.h"
#include "TGraph.h"
#include "TH1.h"
#include "TLine.h"
#include "TMath.h"

using namespace rad;

using std::cout;
using std::endl;
using std::make_unique;
using std::unique_ptr;

/// @brief Calculate the rate of frequency change from radiative emission
/// @param B Magnetic field in Tesla
/// @param theta Pitch angle in radians
/// @param beta Electron speed divided by c
/// @return The rate of change of electron frequency in Hz s^-1
double CalcDfDt(double B, double theta = TMath::PiOver2(),
                double beta = 0.2627) {
  double gamma{1 / sqrt(1 - beta * beta)};
  double E{gamma * ME * TMath::C() * TMath::C()};
  double premult{pow(TMath::Qe(), 5) * pow(B, 3) * TMath::C() /
                 (12 * EPSILON0 * pow(TMath::Pi() * ME, 2))};
  return premult * (1 / (E * E)) * pow(beta * sin(theta), 2) /
         (1 - beta * beta);
}

/// @brief Calculate the optimum acquisition time
/// @param B Magnetic field in Tesla
/// @param theta Pitch angle in radians
/// @param beta Electron speed divided by c
/// @return The optimum acquisition time (for a frequency bin) in seconds
double CalcTAcqOpt(double B, double theta = TMath::PiOver2(),
                   double beta = 0.2627) {
  return 1 / sqrt(CalcDfDt(B, theta, beta));
}

/// @brief Calculate the required power threshold to ensure a given false
/// trigger time
/// @param B Magnetic field in Tesla
/// @param tFalse Mean time (in seconds) between false triggers
/// @param TNoise Noise temperature in kelvin
/// @param obsWidth Observation width in Hz
/// @return The required trigger threshold (in watts)
long double CalcRequiredTriggerThreshold(double B, double TNoise, double tFalse,
                                         double obsWidth) {
  double tAcqOpt{CalcTAcqOpt(B)};
  const double deltaFOpt{1 / tAcqOpt};
  const long double sigma{sqrt(TMath::K() * TNoise * deltaFOpt)};
  const long double nBins{obsWidth / deltaFOpt};
  long double pTrig{tAcqOpt / tFalse};
  long double CDF{pow(1 - pTrig, 1 / nBins)};
  long double vSq{-2 * sigma * sigma * log(1 - CDF)};
  return vSq / 2;
}

/// @brief Calculate the required power threshold to ensure a given false
/// trigger time. Assumes triggering on neighbouring bins
/// @param B Magnetic field in Tesla
/// @param tFalse Mean time (in seconds) between false triggers
/// @param TNoise Noise temperature in kelvin
/// @param obsWidth Observation width in Hz
/// @param nSamples Number of time samples to search in
/// @return The required trigger threshold (in watts)
long double CalcReqTrigThresh(double B, double TNoise, double tFalse,
                              double obsWidth, int nSamples) {
  double tAcqOpt{CalcTAcqOpt(B)};
  const double deltaFOpt{1 / tAcqOpt};
  const long double sigma{sqrt(TMath::K() * TNoise * deltaFOpt)};
  const long double nBins{obsWidth / deltaFOpt};
  long double pTrig{tAcqOpt / tFalse};
  long double CDF{pow(1 - pTrig, 1 / nBins)};
  long double vSq{-2 * sigma * sigma * log(1 - CDF)};
  return vSq / 2;
}

/// @brief Calculates the Larmor power for an electron
/// @param beta Electron speed divided by c
/// @param B Magnetic field in tesla
/// @return Radiated power in watts
double CalcLarmorPower(double beta, double B) {
  double f0{TMath::Qe() * B / ME * (1 / (2 * TMath::Pi()))};
  double premult{2 * TMath::Pi() * pow(TMath::Qe() * f0, 2) /
                 (3 * EPSILON0 * TMath::C())};
  return premult * beta * beta / (1 - beta * beta);
}

unique_ptr<TGraph> MakeTriggerThresholdGraph(double TNoise, double tFalse,
                                             double obsWidth, int nPnts = 200,
                                             double bMin = 0.4,
                                             double bMax = 1.2) {
  auto grThresh = make_unique<TGraph>();
  setGraphAttr(grThresh);
  grThresh->SetLineWidth(3);
  grThresh->SetTitle(Form("T_{noise} = %.0f K; B [T]; P_{trig} [fW]", TNoise));

  for (int n{0}; n < nPnts; n++) {
    double b{bMin + (bMax - bMin) * double(n) / double(nPnts - 1)};
    long double threshold{
        CalcRequiredTriggerThreshold(b, TNoise, tFalse, obsWidth)};
    grThresh->SetPoint(n, b, threshold * 1e15);
  }
  return grThresh;
}

unique_ptr<TGraph> MakeLarmorPowerGraph(int nPnts = 200, double bMin = 0.4,
                                        double bMax = 1.2) {
  auto grLarmor = make_unique<TGraph>();
  setGraphAttr(grLarmor);
  grLarmor->SetTitle("; B [T]; P_{rad} [fW]");
  grLarmor->SetLineWidth(3);

  for (int n{0}; n < nPnts; n++) {
    double b{bMin + (bMax - bMin) * double(n) / double(nPnts - 1)};
    double power{CalcLarmorPower(0.2627, b)};
    grLarmor->SetPoint(n, b, power * 1e15);
  }
  return grLarmor;
}

unique_ptr<TGraph> MakeTriggerThresholdGraphNorm(double TNoise, double tFalse,
                                                 double obsWidth,
                                                 int nPnts = 200,
                                                 double bMin = 0.4,
                                                 double bMax = 1.2) {
  auto grThresh =
      MakeTriggerThresholdGraph(TNoise, tFalse, obsWidth, nPnts, bMin, bMax);
  setGraphAttr(grThresh);
  grThresh->SetLineWidth(3);
  grThresh->SetTitle(
      Form("T_{noise} = %.0f K; B [T]; P_{trig} / P_{rad}", TNoise));

  auto grLarmor = MakeLarmorPowerGraph(nPnts, bMin, bMax);
  for (int n{0}; n < nPnts; n++) {
    double powerNorm{grThresh->GetPointY(n) / grLarmor->GetPointY(n)};
    grThresh->SetPointY(n, powerNorm);
  }
  return grThresh;
}

unique_ptr<TGraph> MakeTrigThreshNeighbour(double TNoise, double tFalse,
                                           double obsWidth, int nSamples,
                                           int nPnts = 200, double bMin = 0.4,
                                           double bMax = 1.2) {
  auto grThresh = make_unique<TGraph>();
  setGraphAttr(grThresh);
  grThresh->SetLineWidth(3);
  grThresh->SetTitle(Form("T_{noise} = %.0f K; B [T]; P_{trig} [fW]", TNoise));

  for (int n{0}; n < nPnts; n++) {
    double b{bMin + (bMax - bMin) * double(n) / double(nPnts - 1)};
    long double threshold{
        CalcRequiredTriggerThreshold(b, TNoise, tFalse, obsWidth)};
    grThresh->SetPoint(n, b, threshold * 1e15);
  }
  return grThresh;
}

unique_ptr<TGraph> MakeTrigThreshNeigbourNorm(double TNoise, double tFalse,
                                              double obsWidth, int nSamples,
                                              int nPnts = 200,
                                              double bMin = 0.4,
                                              double bMax = 1.2) {
  auto grThresh =
      MakeTriggerThresholdGraph(TNoise, tFalse, obsWidth, nPnts, bMin, bMax);
  setGraphAttr(grThresh);
  grThresh->SetLineWidth(3);
  grThresh->SetTitle(
      Form("T_{noise} = %.0f K; B [T]; P_{trig} / P_{rad}", TNoise));

  auto grLarmor = MakeLarmorPowerGraph(nPnts, bMin, bMax);
  for (int n{0}; n < nPnts; n++) {
    double powerNorm{grThresh->GetPointY(n) / grLarmor->GetPointY(n)};
    grThresh->SetPointY(n, powerNorm);
  }
  return grThresh;
}

unique_ptr<TGraph> MakeTfPlot(double T, double bandwidth, int nBins = 1,
                              int nPnts = 400, double pMin = 5e-18,
                              double pMax = 0.5e-15) {
  const long double sigma{sqrt(TMath::K() * T * bandwidth)};
  const long double tAcq{1 / bandwidth};
  const double logPowerDiff{(log10(pMax) - log10(pMin)) / double(nPnts - 1)};

  auto grTime = make_unique<TGraph>();
  setGraphAttr(grTime);
  const long double nBins100MHz{100e6 / bandwidth};
  grTime->SetLineWidth(3);
  grTime->SetTitle(Form(
      "T_{noise} = %.0f K; P_{bin} [fW]; Time between false triggers [s]", T));

  for (int n{0}; n < nPnts; n++) {
    double power{pMin * pow(10, logPowerDiff * double(n))};
    long double vMag{sqrt(2 * power)};
    long double probNoise100MHz{1 - pow(RayleighCDF(vMag, sigma), nBins100MHz)};
    long double pTrig{pow(probNoise100MHz, nBins)};
    if (pTrig > 0) {
      long double trigTime{tAcq / pTrig};
      grTime->SetPoint(n, power * 1e15, trigTime);
    }
  }
  return grTime;
}

unique_ptr<TGraph> MakeTfPlotNeighbour(double T, double bandwidth, int nBins,
                                       int nPnts = 400, double pMin = 5e-18,
                                       double pMax = 0.5e-15) {
  const long double sigma{sqrt(TMath::K() * T * bandwidth)};
  const long double tAcq{1 / bandwidth};
  const double logPowerDiff{(log10(pMax) - log10(pMin)) / double(nPnts - 1)};

  auto grTime = make_unique<TGraph>();
  setGraphAttr(grTime);
  const long double nBins100MHz{100e6 / bandwidth};
  grTime->SetLineWidth(3);
  grTime->SetTitle(Form(
      "T_{noise} = %.0f K; P_{bin} [fW]; Time between false triggers [s]", T));

  for (int n{0}; n < nPnts; n++) {
    double power{pMin * pow(10, logPowerDiff * double(n))};
    long double vMag{sqrt(2 * power)};
    long double probNoise100MHz{1 - pow(RayleighCDF(vMag, sigma), nBins100MHz)};
    long double probNoise1Bin{1 - RayleighCDF(vMag, sigma)};
    long double pTrig{probNoise100MHz * pow(probNoise1Bin, nBins - 1) *
                      pow(2, nBins - 1)};
    if (pTrig > 0) {
      long double trigTime{tAcq / pTrig};
      grTime->SetPoint(n, power * 1e15, trigTime);
    }
  }
  return grTime;
}

unique_ptr<TGraph> MakeFalseTrigProbPlot(double T, double bandwidth,
                                         int nPnts = 400, double pMin = 5e-18,
                                         double pMax = 0.1e-15) {
  const long double sigma{sqrt(TMath::K() * T * bandwidth)};
  const double logPowerDiff{(log10(pMax) - log10(pMin)) / double(nPnts - 1)};

  auto grProb = make_unique<TGraph>();
  setGraphAttr(grProb);
  const long double nBins100MHz{100e6 / bandwidth};
  grProb->SetLineWidth(3);
  grProb->SetTitle(Form("T_{noise} = %.0f K; P_{bin} [fW]; P (false trig)", T));

  for (int n{0}; n < nPnts; n++) {
    double power{pMin * pow(10, logPowerDiff * double(n))};
    long double vMag{sqrt(2 * power)};
    long double probNoise100MHz{1 - pow(RayleighCDF(vMag, sigma), nBins100MHz)};
    if (probNoise100MHz > 0) {
      grProb->SetPoint(n, power * 1e15, probNoise100MHz);
    }
  }
  return grProb;
}

int main(int argc, char *argv[]) {
  TString outputFile{argv[1]};
  auto fout{make_unique<TFile>(outputFile, "RECREATE")};

  // Set up noise parameters
  const double noiseTemp{11};    // Kelvin
  const double bandwidth{17e3};  // Hz
  const long double sigma{sqrt(TMath::K() * noiseTemp * bandwidth)};
  cout << "Sigma = " << sigma << endl;
  const long double mean{sigma * sqrt(TMath::Pi() / 2)};
  cout << "Power per bandwidth = " << TMath::K() * noiseTemp << " W/Hz\n";
  cout << "Power per 17 kHz bin = " << TMath::K() * noiseTemp * bandwidth
       << " W\n";

  const long double maxDistX{6 * sigma};
  const int nFuncPnts{10000};
  auto fCDF{make_unique<TF1>("fCDF", RayleighCDFFunc, 0, maxDistX, 1)};
  fCDF->SetParameter(0, sigma);
  fCDF->SetNpx(nFuncPnts);
  fCDF->SetLineWidth(3);
  fCDF->SetTitle(Form("CDF, #sigma = %.3Le; x; F(x; #sigma)", sigma));
  auto fPDF{make_unique<TF1>("fPDF", RayleighPDFFunc, 0, maxDistX, 1)};
  fPDF->SetParameter(0, sigma);
  fPDF->SetNpx(nFuncPnts);
  fPDF->SetLineWidth(3);
  fPDF->SetTitle(Form("PDF, #sigma = %.3Le; x; f(x; #sigma)", sigma));

  const double maxXExample{10};
  std::vector<double> sigmaExamples{0.5, 1, 2, 3, 4};
  for (auto &s : sigmaExamples) {
    auto fPDFEx = make_unique<TF1>(Form("#sigma = %.1f", s), RayleighPDFFunc, 0,
                                   maxXExample, 1);
    fPDFEx->SetParameter(0, s);
    fPDFEx->SetLineWidth(3);
    fPDFEx->GetXaxis()->SetTitle("x");
    fPDFEx->GetYaxis()->SetTitle("f(x; #sigma)");
    auto fCDFEx = make_unique<TF1>(Form("#sigma = %.1f", s), RayleighCDFFunc, 0,
                                   maxXExample, 1);
    fCDFEx->SetParameter(0, s);
    fCDFEx->SetLineWidth(3);
    fCDFEx->GetXaxis()->SetTitle("x");
    fCDFEx->GetYaxis()->SetTitle("F(x; #sigma)");

    fout->cd();
    fPDFEx->Write(Form("fPDFEx_%.0f", s));
    fCDFEx->Write(Form("fCDFEx_%.0f", s));
  }

  auto meanLine{make_unique<TLine>(mean, 0, mean, fPDF->Eval(mean))};
  meanLine->SetLineWidth(2);
  meanLine->SetLineStyle(2);

  fout->cd();
  fCDF->Write("fCDF");
  fPDF->Write("fPDF");
  meanLine->Write("meanLine");

  const int nBins{5882};
  auto hWhiteNoise =
      make_unique<TH1D>("hWhiteNoise", "; Frequency [MHz]; FFT [W^{0.5}]",
                        nBins, 0, nBins * bandwidth / 1e6);
  SetHistAttr(hWhiteNoise);
  auto hFFT =
      make_unique<TH1D>("hFFT", "; FFT [W^{0.5}]; N_{bins}", 200, 0, 6 * sigma);
  SetHistAttr(hFFT);

  double totalPower{0};
  for (int n{0}; n < nBins; n++) {
    double fft{fPDF->GetRandom()};
    hFFT->Fill(fft);
    hWhiteNoise->SetBinContent(n + 1, fft);
    double power{fft * fft / 2};
    totalPower += power;
  }
  double sampledPowerPerBin{totalPower / double(nBins)};
  cout << "Sampled power per bin = " << sampledPowerPerBin << " W\n";

  fout->cd();
  hFFT->Write();
  hWhiteNoise->Write();
  hFFT.reset();
  hWhiteNoise.reset();

  const int nBins100MHz{5882};
  const double hypotheticalSignalPower{1e-16};  // W
  const long double sigVoltageMag{sqrt(2 * hypotheticalSignalPower)};
  const long double probNoiseHigher1Bin{1 - RayleighCDF(sigVoltageMag, sigma)};
  const long double probNoiseHigher100MHz{
      1 - pow(RayleighCDF(sigVoltageMag, sigma), nBins100MHz)};
  cout << "Probability a given noise bin is higher = " << probNoiseHigher1Bin
       << endl;
  cout << "Probability that any noise bin in a 100 MHz "
          "range is higher = "
       << probNoiseHigher100MHz << endl;

  // Plot probabilities as a function of bin power
  // Also do this for different noise temperatures
  auto grProb6K = MakeFalseTrigProbPlot(6, bandwidth);
  grProb6K->SetLineColor(kRed);
  auto grProb7K = MakeFalseTrigProbPlot(7, bandwidth);
  grProb7K->SetLineColor(kBlue);
  auto grProb11K = MakeFalseTrigProbPlot(11, bandwidth);
  grProb11K->SetLineColor(kOrange + 1);

  auto grTime6K = MakeTfPlot(6, bandwidth, 1);
  grTime6K->SetLineColor(kRed);
  auto grTime7K = MakeTfPlot(7, bandwidth, 1);
  grTime7K->SetLineColor(kBlue);
  auto grTime11K = MakeTfPlot(11, bandwidth, 1);
  grTime11K->SetLineColor(kOrange + 1);

  auto grTime6K_2bins = MakeTfPlot(6, bandwidth, 2);
  grTime6K_2bins->SetLineColor(kMagenta + 1);
  auto grTime7K_2bins = MakeTfPlot(7, bandwidth, 2);
  grTime7K_2bins->SetLineColor(kMagenta + 1);
  auto grTime11K_2bins = MakeTfPlot(11, bandwidth, 2);
  grTime11K_2bins->SetLineColor(kMagenta + 1);

  auto grTime6K_4bins = MakeTfPlot(6, bandwidth, 4);
  grTime6K_4bins->SetLineColor(kBlack);
  auto grTime7K_4bins = MakeTfPlot(7, bandwidth, 4);
  grTime7K_4bins->SetLineColor(kBlack);
  auto grTime11K_4bins = MakeTfPlot(11, bandwidth, 4);
  grTime11K_4bins->SetLineColor(kBlack);

  auto grTime6K_2binsC = MakeTfPlotNeighbour(6, bandwidth, 2);
  grTime6K_2binsC->SetLineColor(kGreen + 1);
  auto grTime7K_2binsC = MakeTfPlotNeighbour(7, bandwidth, 2);
  grTime7K_2binsC->SetLineColor(kGreen + 1);
  auto grTime11K_2binsC = MakeTfPlotNeighbour(11, bandwidth, 2);
  grTime11K_2binsC->SetLineColor(kGreen + 1);

  auto grTime6K_4binsC = MakeTfPlotNeighbour(6, bandwidth, 4);
  grTime6K_4binsC->SetLineColor(kCyan + 1);
  auto grTime7K_4binsC = MakeTfPlotNeighbour(7, bandwidth, 4);
  grTime7K_4binsC->SetLineColor(kCyan + 1);
  auto grTime11K_4binsC = MakeTfPlotNeighbour(11, bandwidth, 4);
  grTime11K_4binsC->SetLineColor(kCyan + 1);

  // 0.7 T
  auto grTime11K_700mT_4binsC = MakeTfPlotNeighbour(11, 11e3, 4);
  grTime11K_700mT_4binsC->SetLineColor(kCyan + 1);

  auto grLarmor = MakeLarmorPowerGraph();
  grLarmor->SetLineColor(kMagenta + 1);
  grLarmor->SetLineStyle(2);

  auto grThresh6K = MakeTriggerThresholdGraph(6, 10, 100e6);
  grThresh6K->SetLineColor(kRed);
  auto grThresh7K = MakeTriggerThresholdGraph(7, 10, 100e6);
  grThresh7K->SetLineColor(kBlue);
  auto grThresh11K = MakeTriggerThresholdGraph(11, 10, 100e6);
  grThresh11K->SetLineColor(kOrange + 1);

  auto grThreshNorm6K = MakeTriggerThresholdGraphNorm(6, 10, 100e6);
  grThreshNorm6K->SetLineColor(kRed);
  auto grThreshNorm7K = MakeTriggerThresholdGraphNorm(7, 10, 100e6);
  grThreshNorm7K->SetLineColor(kBlue);
  auto grThreshNorm11K = MakeTriggerThresholdGraphNorm(11, 10, 100e6);
  grThreshNorm11K->SetLineColor(kOrange + 1);

  /*
  auto grThreshNorm7K_4binsC = MakeTrigThreshNeighbourNorm(6, 10, 100e6, 4);
  grThreshNorm7K_4binsC->SetLineColor(kMagenta + 1);
  auto grThreshNorm7K_4binsC = MakeTrigThreshNeighbourNorm(6, 10, 100e6, 4);
  grThreshNorm7K_4binsC->SetLineColor(kMagenta + 1);
  */

  fout->cd();
  grProb6K->Write("grProb6K");
  grProb7K->Write("grProb7K");
  grProb11K->Write("grProb11K");

  grTime6K->Write("grTime6K");
  grTime7K->Write("grTime7K");
  grTime11K->Write("grTime11K");

  grTime6K_2bins->Write("grTime6K_2bins");
  grTime7K_2bins->Write("grTime7K_2bins");
  grTime11K_2bins->Write("grTime11K_2bins");

  grTime6K_4bins->Write("grTime6K_4bins");
  grTime7K_4bins->Write("grTime7K_4bins");
  grTime11K_4bins->Write("grTime11K_4bins");

  grTime6K_2binsC->Write("grTime6K_2binsC");
  grTime7K_2binsC->Write("grTime7K_2binsC");
  grTime11K_2binsC->Write("grTime11K_2binsC");

  grTime6K_4binsC->Write("grTime6K_4binsC");
  grTime7K_4binsC->Write("grTime7K_4binsC");
  grTime11K_4binsC->Write("grTime11K_4binsC");
  grTime11K_700mT_4binsC->Write("grTime11K_700mT_4binsC");

  grLarmor->Write("grLarmor");
  grThresh6K->Write("grThresh6K");
  grThresh7K->Write("grThresh7K");
  grThresh11K->Write("grThresh11K");

  grThreshNorm6K->Write("grThreshNorm6K");
  grThreshNorm7K->Write("grThreshNorm7K");
  grThreshNorm11K->Write("grThreshNorm11K");

  fout->Close();
  return 0;
}