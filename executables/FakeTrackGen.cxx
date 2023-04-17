/*
  FakeTrackGen.cxx

  Make some fake tracks with a given SNR
*/

#include <chrono>
#include <random>

#include "BasicFunctions/BasicFunctions.h"
#include "BasicFunctions/Constants.h"
#include "TFile.h"
#include "TGraph.h"
#include "TH2.h"
#include "TMath.h"
#include "TString.h"

using namespace rad;

TGraph *GetChirpGraph(double srate, int N, double A, double f0, double c) {
  TGraph *gr = new TGraph();
  setGraphAttr(gr);
  gr->GetXaxis()->SetTitle("Time [s]");
  gr->GetYaxis()->SetTitle("Amplitude [A. U.]");

  double timeInterval{1.0 / srate};

  // Generate a random starting phase
  long int seed{std::chrono::system_clock::now().time_since_epoch().count()};
  std::mt19937 rng(seed);
  std::uniform_real_distribution<double> dist(0, 2 * TMath::Pi());
  double phi0{dist(rng)};

  for (int iS{0}; iS < N; iS++) {
    double time{double(iS) * timeInterval};
    double x{ChirpSignal(A, time, phi0, f0, c)};
    gr->SetPoint(iS, time, x);
  }

  return gr;
}

void AddWhiteNoise(TGraph *gr, double T) {
  const double fs{1.0 / (gr->GetPointX(1) - gr->GetPointX(0))};
  const double sigma{sqrt(TMath::K() * T * fs / 2)};
  long int seed{std::chrono::system_clock::now().time_since_epoch().count()};
  std::mt19937 rng(seed);
  std::normal_distribution<double> noise{0, sigma};
  for (int n{0}; n < gr->GetN(); n++) {
    double amp{gr->GetPointY(n)};
    amp += noise(rng);
    gr->SetPointY(n, amp);
  }
}

double GetChirpRate(double Ek, double B, double theta) {
  const double ETot{Ek * TMath::Qe() + ME * pow(TMath::C(), 2)};
  const double gamma{ETot / (ME * pow(TMath::C(), 2))};
  const double betaSq{1 - 1 / (gamma * gamma)};
  return betaSq * pow(sin(theta), 2) * pow(TMath::Qe(), 5) * pow(B, 3) *
         TMath::C() /
         (ETot * ETot * (1 - betaSq) * 12 * pow(TMath::Pi() * ME, 2) *
          EPSILON0);
}

TH2D *MakeSpectrogram(double srate, double totTime, double snrMax, double B = 1,
                      double theta = TMath::PiOver2(), double Ek = 18.6e3) {
  // Calculate chirp rate
  const double cRate{GetChirpRate(Ek, B, theta)};

  // Optimal time and frequency bins
  const double tAcqOpt{pow(cRate, -0.5)};
  const double deltaFOpt{1 / tAcqOpt};
  // Calculate the number of time bins
  const int nTimeBins{int(std::floor(totTime / tAcqOpt))};
  // Number of samples per time bin
  const int nSamplesPerTimeBin{int(std::round(tAcqOpt * srate))};
  const double fBinWidth{srate / double(nSamplesPerTimeBin)};
  const int nFBins{int(std::round((srate / 2) / fBinWidth))};

  // Initialize histogram
  auto h2Spec = new TH2D(
      "",
      Form("SNR_{max} = %.0f; Time [ms]; Frequency [MHz]; Power [A. U.]",
           snrMax),
      nTimeBins, 0, double(nTimeBins) * tAcqOpt * 1e3, nFBins, 0,
      double(nFBins) * fBinWidth / 1e6);
  SetHistAttr(h2Spec);

  // Calculate start frequency
  const double f0{CalcCyclotronFreq(Ek, B)};  // Hz
  // Now after downmixing
  const double mixerFreq{26.81e9};    // Hz
  const double f0DM{f0 - mixerFreq};  // Hz
  std::cout << "Downmixed frequency = " << f0DM / 1e6 << " MHz\n";

  // Make the signal graph
  const double noisePower{1};
  const double signalPower{noisePower * snrMax * sin(theta) * sin(theta)};
  auto gr = GetChirpGraph(srate, nSamplesPerTimeBin * nTimeBins,
                          sqrt(2 * signalPower), f0DM, cRate);
  // Add noise
  const double noiseTemp{noisePower / (TMath::K() * deltaFOpt)};
  AddWhiteNoise(gr, noiseTemp);

  // Now divide the graph up and add to the spectrogram
  for (int iBin{1}; iBin <= h2Spec->GetNbinsX(); iBin++) {
    // Create a graph which is the relevant subset of points
    const int p1{(iBin - 1) * nSamplesPerTimeBin};
    const int p2{iBin * nSamplesPerTimeBin};

    auto grSub = new TGraph();
    for (int n{p1}; n < p2; n++) {
      grSub->SetPoint(grSub->GetN(), gr->GetPointX(n), gr->GetPointY(n));
    }
    // Do the FFT
    auto grSubPower = MakePowerSpectrumPeriodogram(grSub);
    delete grSub;

    // Now write to the histogram
    for (int iF{0}; iF < grSubPower->GetN(); iF++) {
      h2Spec->SetBinContent(iBin, iF + 1, grSubPower->GetPointY(iF));
    }
    delete grSubPower;
  }

  return h2Spec;
}

int main(int argc, char *argv[]) {
  TString outputFile{argv[1]};
  auto fout = std::make_unique<TFile>(outputFile, "recreate");

  const double sampleRate{1e9};   // Samples per second
  const double chosenField{1.0};  // Tesla

  long int seed{std::chrono::system_clock::now().time_since_epoch().count()};
  std::mt19937 rng(seed);

  const double runTime{1e-3};  //  seconds

  for (int SNR{1}; SNR <= 18; SNR++) {
    // Create 10 spectograms with slight variations in theta and Ek
    const int nThrows{10};
    std::uniform_real_distribution<double> thetaDegDist(89.5, 90);
    std::uniform_real_distribution<double> EkDist(18.5e3, 18.6e3);
    for (int iTh{0}; iTh < nThrows; iTh++) {
      double thetaDeg{thetaDegDist(rng)};
      double theta{thetaDeg * TMath::Pi() / 180};
      double Ek{EkDist(rng)};
      auto h2Spec = MakeSpectrogram(sampleRate, runTime, double(SNR),
                                    chosenField, theta, Ek);
      h2Spec->SetTitle(
          Form("SNR_{max} = %d, #theta = %.2f^{#circ}, E = %.0f eV", SNR,
               thetaDeg, Ek));

      // Calculate the downmixed frequency
      // Set the axis to the region of interest
      const double f0{CalcCyclotronFreq(Ek, chosenField)};
      const double f0DM{f0 - 26.81e9};
      h2Spec->GetYaxis()->SetRangeUser(f0DM / 1e6 - 1, f0DM / 1e6 + 1);
      fout->cd();
      h2Spec->Write(Form("h2SpecSNR%d_%d", SNR, iTh));
    }
  }

  fout->Close();
  return 0;
}