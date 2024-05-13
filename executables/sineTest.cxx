// sineTest.cxx

// STL
#include <cmath>
#include <iostream>

// ROOT includes
#include "BasicFunctions/BasicFunctions.h"
#include "FieldClasses/FieldClasses.h"
#include "TAxis.h"
#include "TFile.h"
#include "TGraph.h"
#include "TMath.h"
#include "TRandom3.h"

int main() {
  double totalTime1 = 500.0;
  double stepSize1 = 0.01;
  int nSteps1 = totalTime1 / stepSize1;
  double totalTime2 = 300.0;
  double stepSize2 = 0.005;
  int nSteps2 = totalTime2 / stepSize2;

  TRandom3* randomNum = new TRandom3();

  TGraph* gr1 = new TGraph();
  for (int n = 0; n < nSteps1; n++) {
    double time = n * stepSize1;
    double amp = TMath::Sin(2 * TMath::Pi() * 1.0 * time) +
                 0.2 * TMath::Sin(2 * TMath::Pi() * 5.0 * time + 1.0) +
                 0.3 * TMath::Sin(2 * TMath::Pi() * 7.0 * time + 1.5) +
                 randomNum->Gaus(0, 0.05);
    gr1->SetPoint(gr1->GetN(), time, amp);
  }

  TGraph* gr2 = new TGraph();
  for (int n = 0; n < nSteps2; n++) {
    double time = n * stepSize2;
    double amp = TMath::Sin(2 * TMath::Pi() * 1.0 * time) +
                 0.2 * TMath::Sin(2 * TMath::Pi() * 5.0 * time + 1.0) +
                 0.3 * TMath::Sin(2 * TMath::Pi() * 7.0 * time + 1.5) +
                 randomNum->Gaus(0, 0.05);
    gr2->SetPoint(gr2->GetN(), time, amp);
  }

  TGraph* gr1FFT = rad::MakePowerSpectrumPeriodogram(gr1);
  TGraph* gr2FFT = rad::MakePowerSpectrumPeriodogram(gr2);
  TGraph* gr1FFTVs = rad::MakePowerSpectrumNorm(gr1);
  TGraph* gr2FFTVs = rad::MakePowerSpectrumNorm(gr2);
  gr1FFT->SetLineWidth(2);
  gr2FFT->SetLineWidth(2);
  gr1FFTVs->SetLineWidth(2);
  gr2FFTVs->SetLineWidth(2);

  std::cout << "FFT integral (1, 2) = " << rad::IntegratePowerNorm(gr1FFTVs)
            << ", " << rad::IntegratePowerNorm(gr2FFTVs) << std::endl;

  TFile* fout = new TFile("outputSine.root", "recreate");
  gr1->Write("gr1");
  gr2->Write("gr2");
  gr1FFT->Write("gr1FFT");
  gr2FFT->Write("gr2FFT");
  gr1FFTVs->Write("gr1FFTVs");
  gr2FFTVs->Write("gr2FFTVs");

  fout->Close();
  delete fout;

  return 0;
}
