/*
  SignalRequirements.cxx

  Using known noise properties, calculate the required signal power
  Allow for some hypothesised noise rate
*/

#include <iostream>
#include <memory>
#include <span>

#include "BasicFunctions/BasicFunctions.h"
#include "BasicFunctions/Constants.h"
#include "TF1.h"
#include "TFile.h"
#include "TGraph.h"
#include "TLine.h"
#include "TMath.h"

using namespace rad;

using std::cout;
using std::endl;
using std::make_unique;
using std::unique_ptr;

unique_ptr<TGraph> MakeFalseTrigTimePlot(double T, double bandwidth,
                                         int nPnts = 400, double pMin = 5e-18,
                                         double pMax = 1e-15) {
  const long double sigma{TMath::K() * T * bandwidth};
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
    if (probNoise100MHz > 0) {
      long double trigTime{tAcq / probNoise100MHz};
      grTime->SetPoint(n, power * 1e15, trigTime);
    }
  }
  return grTime;
}

unique_ptr<TGraph> MakeFalseTrigProbPlot(double T, double bandwidth,
                                         int nPnts = 400, double pMin = 5e-18,
                                         double pMax = 1e-15) {
  const long double sigma{TMath::K() * T * bandwidth};
  const long double tAcq{1 / bandwidth};
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
  auto args = std::span(argv, size_t(argc));
  TString outputFile{args[1]};
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

  const int nBins{50000};
  double totalPower{0};
  for (int n{0}; n < nBins; n++) {
    double power{pow(fPDF->GetRandom(), 2) / 2};
    totalPower += power;
  }
  double sampledPowerPerBin{totalPower / double(nBins)};
  cout << "Sampled power per bin = " << sampledPowerPerBin << " W\n";

  const int nBins100MHz{5882};
  const double hypotheticalSignalPower{1e-16};  // W
  const long double sigVoltageMag{sqrt(2 * hypotheticalSignalPower)};
  const long double probNoiseHigher1Bin{1 - RayleighCDF(sigVoltageMag, sigma)};
  const long double probNoiseHigher100MHz{
      1 - pow(RayleighCDF(sigVoltageMag, sigma), nBins100MHz)};
  cout << "Probability a given noise bin is higher = " << probNoiseHigher1Bin
       << endl;
  cout << "Probability that any noise bin in a 100 MHz range is higher = "
       << probNoiseHigher100MHz << endl;

  // Plot probabilities as a function of bin power
  // Also do this for different noise temperatures
  auto grProb6K = MakeFalseTrigProbPlot(6, bandwidth);
  grProb6K->SetLineColor(kRed);
  auto grProb7K = MakeFalseTrigProbPlot(7, bandwidth);
  grProb7K->SetLineColor(kBlue);
  auto grProb11K = MakeFalseTrigProbPlot(11, bandwidth);
  grProb11K->SetLineColor(kOrange + 1);

  auto grTime6K = MakeFalseTrigTimePlot(6, bandwidth);
  grTime6K->SetLineColor(kRed);
  auto grTime7K = MakeFalseTrigTimePlot(7, bandwidth);
  grTime7K->SetLineColor(kBlue);
  auto grTime11K = MakeFalseTrigTimePlot(11, bandwidth);
  grTime11K->SetLineColor(kOrange + 1);

  fout->cd();
  grProb6K->Write("grProb6K");
  grProb7K->Write("grProb7K");
  grProb11K->Write("grProb11K");
  grTime6K->Write("grTime6K");
  grTime7K->Write("grTime7K");
  grTime11K->Write("grTime11K");

  fout->Close();
  return 0;
}