/// ScaledSignal.cxx

#include "physics/SignalProcessing/ScaledSignal.h"

#include <cassert>
#include <cmath>

rad::ScaledSignal::ScaledSignal(InducedVoltage iv, LocalOscillator lo,
                                double srate, double scaleFac,
                                std::vector<GaussianNoise> noiseTerms,
                                double maxTime) {
  assert(!isnan(scaleFac) && !isinf(scaleFac));
  sampleRate = srate;
  scaleFactor = scaleFac;

  // Make sure the noise terms are all set up correctly
  for (int iNoise = 0; iNoise < noiseTerms.size(); iNoise++) {
    (noiseTerms.at(iNoise)).SetSampleFreq(sampleRate);
    (noiseTerms.at(iNoise)).SetSigma();
  }

  // The actual output graphs
  grVITime = new TGraph();
  grVQTime = new TGraph();
  setGraphAttr(grVITime);
  setGraphAttr(grVQTime);
  grVITime->GetYaxis()->SetTitle("V_{I}");
  grVQTime->GetYaxis()->SetTitle("V_{Q}");
  grVITime->GetXaxis()->SetTitle("Time [s]");
  grVQTime->GetXaxis()->SetTitle("Time [s]");

  // Split the signal up into chunks to avoid memory issues
  if (maxTime < 0) maxTime = iv.GetFinalTime();

  const double chunkSize = 25e-6;
  double lastChunk = 0;
  double thisChunk = lastChunk + chunkSize;
  if (thisChunk > maxTime) thisChunk = maxTime;

  double thisSample = 0;
  double this10Sample = 0;

  while (thisChunk <= maxTime && thisChunk != lastChunk) {
    ProcessTimeChunk(iv, lo, thisChunk, lastChunk, noiseTerms, thisSample,
                     this10Sample);
    lastChunk = thisChunk;
    thisChunk += chunkSize;
    if (thisChunk > maxTime) thisChunk = maxTime;
  }
}

void rad::ScaledSignal::ProcessTimeChunk(InducedVoltage iv, LocalOscillator lo,
                                         double thisChunk, double lastChunk,
                                         std::vector<GaussianNoise> noiseTerms,
                                         double& firstSampleTime,
                                         double& firstSample10Time) {
  iv.ResetVoltage();
  iv.GenerateVoltage(lastChunk, thisChunk);

  TGraph* grInputVoltageTemp = iv.GetVoltageGraph();

  // Scale the input voltage
  ScaleGraph(grInputVoltageTemp, scaleFactor);

  std::cout << "Performing the downmixing..." << std::endl;
  TGraph* grVITimeUnfiltered = DownmixInPhase(grInputVoltageTemp, lo);
  TGraph* grVQTimeUnfiltered = DownmixQuadrature(grInputVoltageTemp, lo);

  std::cout << "First downsampling" << std::endl;
  TGraph* grVITimeFirstSample =
      SampleWaveform(grVITimeUnfiltered, 10 * sampleRate, firstSample10Time);
  TGraph* grVQTimeFirstSample =
      SampleWaveform(grVQTimeUnfiltered, 10 * sampleRate, firstSample10Time);
  delete grVITimeUnfiltered;
  delete grVQTimeUnfiltered;
  firstSample10Time =
      grVITimeFirstSample->GetPointX(grVITimeFirstSample->GetN() - 1) +
      1 / (10 * sampleRate);

  // Now we need to filter and then sample these signals
  std::cout << "Filtering.." << std::endl;
  TGraph* grVITimeUnsampled =
      BandPassFilter(grVITimeFirstSample, 0.0, sampleRate / 2.0);
  TGraph* grVQTimeUnsampled =
      BandPassFilter(grVQTimeFirstSample, 0.0, sampleRate / 2.0);

  delete grVITimeFirstSample;
  delete grVQTimeFirstSample;

  // Now do sampling
  // Use simple linear interpolation for the job
  std::cout << "Sampling..." << std::endl;
  TGraph* grVITimeTemp =
      SampleWaveform(grVITimeUnsampled, sampleRate, firstSampleTime);
  TGraph* grVQTimeTemp =
      SampleWaveform(grVQTimeUnsampled, sampleRate, firstSampleTime);
  delete grVITimeUnsampled;
  delete grVQTimeUnsampled;
  firstSampleTime =
      grVITimeTemp->GetPointX(grVITimeTemp->GetN() - 1) + 1 / sampleRate;

  std::cout << "Adding noise..." << std::endl;
  AddGaussianNoise(grVITimeTemp, noiseTerms);
  AddGaussianNoise(grVQTimeTemp, noiseTerms);

  if (iv.GetLowerAntennaBandwidth() != -DBL_MAX ||
      iv.GetUpperAntennaBandwidth() != DBL_MAX) {
    std::cout << "Implementing antenna bandwidth..." << std::endl;
    grVITimeTemp = BandPassFilter(
        grVITimeTemp, iv.GetLowerAntennaBandwidth() - lo.GetFrequency(),
        iv.GetUpperAntennaBandwidth() - lo.GetFrequency());
    grVQTimeTemp = BandPassFilter(
        grVQTimeTemp, iv.GetLowerAntennaBandwidth() - lo.GetFrequency(),
        iv.GetUpperAntennaBandwidth() - lo.GetFrequency());
  }

  // Now add the information from these temporary graphs to the larger ones
  for (int i = 0; i < grVITimeTemp->GetN(); i++) {
    grVITime->SetPoint(grVITime->GetN(), grVITimeTemp->GetPointX(i),
                       grVITimeTemp->GetPointY(i));
    grVQTime->SetPoint(grVQTime->GetN(), grVQTimeTemp->GetPointX(i),
                       grVQTimeTemp->GetPointY(i));
  }

  delete grVITimeTemp;
  delete grVQTimeTemp;
}

rad::ScaledSignal::ScaledSignal(const ScaledSignal& s1) {
  sampleRate = s1.sampleRate;
  scaleFactor = s1.scaleFactor;
  grVITime = (TGraph*)s1.grVITime->Clone();
  grVQTime = (TGraph*)s1.grVQTime->Clone();
}
