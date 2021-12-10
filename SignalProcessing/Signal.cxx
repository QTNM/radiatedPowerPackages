// Signal.cxx

#include "SignalProcessing/Signal.h"
#include "SignalProcessing/InducedVoltage.h"
#include "SignalProcessing/LocalOscillator.h"
#include "FieldClasses/FieldClasses.h"
#include "BasicFunctions/BasicFunctions.h"

#include <vector>
#include <iostream>
#include <cassert>
#include <cmath>

#include "TGraph.h"
#include "TH2.h"

#include "FFTtools.h"

rad::Signal::~Signal() {
  delete grVITime;
  delete grVQTime;
}

rad::Signal::Signal(std::vector<FieldPoint> fp, LocalOscillator lo, double srate,
		    std::vector<GaussianNoise> noiseTerms, const bool kUseRetardedTime) {
  assert(fp.size() > 0);
  
  sampleRate = srate;
  
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
  grVITime->GetXaxis()->SetTitle("Time [s]");
  grVQTime->GetXaxis()->SetTitle("Time [s]");
  
  std::vector<TGraph*> vecVITime;
  std::vector<TGraph*> vecVQTime;
  
  // For each antenna point generate the field and do the signal processing
  for (int point = 0; point < fp.size(); point++) {
    TGraph* grInputVoltageTemp = fp[point].GetAntennaLoadVoltageTimeDomain(kUseRetardedTime, -1, -1);
    
    std::cout<<"Performing the downmixing..."<<std::endl;
    TGraph* grVITimeUnfiltered = DownmixInPhase(grInputVoltageTemp, lo);
    TGraph* grVQTimeUnfiltered = DownmixQuadrature(grInputVoltageTemp, lo);

    std::cout<<"First downsampling..."<<std::endl;
    TGraph* grVITimeFirstSample = SampleWaveform(grVITimeUnfiltered, 10*sampleRate);
    TGraph* grVQTimeFirstSample = SampleWaveform(grVQTimeUnfiltered, 10*sampleRate);

    delete grVITimeUnfiltered;
    delete grVQTimeUnfiltered;
    
    // Now we need to filter and then sample these signals
    std::cout<<"Filtering.."<<std::endl;
    TGraph* grVITimeUnsampled = BandPassFilter(grVITimeFirstSample, 0.0, sampleRate/2.0);
    TGraph* grVQTimeUnsampled = BandPassFilter(grVQTimeFirstSample, 0.0, sampleRate/2.0);

    delete grVITimeFirstSample;
    delete grVQTimeFirstSample;
    
    // Now do sampling
    // Use simple linear interpolation for the job
    std::cout<<"Sampling..."<<std::endl;
    TGraph* grVITimeTemp = SampleWaveform(grVITimeUnsampled);
    TGraph* grVQTimeTemp = SampleWaveform(grVQTimeUnsampled);

    delete grVITimeUnsampled;
    delete grVQTimeUnsampled;
    
    std::cout<<"Adding noise..."<<std::endl;
    AddGaussianNoise(grVITimeTemp, noiseTerms);
    AddGaussianNoise(grVQTimeTemp, noiseTerms);

    vecVITime.push_back(grVITimeTemp);
    vecVQTime.push_back(grVQTimeTemp);
  } // Generate individual waveforms

  std::cout<<"Summing fields"<<std::endl;
  const int nSampledTimePoints = vecVITime[0]->GetN();
  for (int t = 0; t < nSampledTimePoints; t++) {
    // Loop over the signals from each individual antenna and sum
    double sumVI = 0;
    double sumVQ = 0;
    for (int field = 0; field < vecVITime.size(); field++) {
      sumVI += vecVITime[field]->GetPointY(t);
      sumVQ += vecVQTime[field]->GetPointY(t);     
    }
    grVITime->SetPoint(t, vecVITime[0]->GetPointX(t), sumVI);
    grVQTime->SetPoint(t, vecVQTime[0]->GetPointX(t), sumVQ);
  }
  
  // Remaining cleanup
  // for (auto p : vecVITime) {
  //    delete p;
  // }
  // vecVITime.clear();
  
  // for (auto p : vecVQTime) {
  //    delete p;
  // }
  // vecVQTime.clear();
  
  std::cout<<"Constructor completed"<<std::endl;
  
}

rad::Signal::Signal(FieldPoint fp, LocalOscillator lo, double srate,
		    std::vector<GaussianNoise> noiseTerms, const bool kUseRetardedTime) {
  sampleRate = srate;
  
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
  grVITime->GetXaxis()->SetTitle("Time [s]");
  grVQTime->GetXaxis()->SetTitle("Time [s]");

  TGraph* grInputVoltageTemp = fp.GetAntennaLoadVoltageTimeDomain(kUseRetardedTime, -1, -1);
    
  std::cout<<"Performing the downmixing..."<<std::endl;
  TGraph* grVITimeUnfiltered = DownmixInPhase(grInputVoltageTemp, lo);
  TGraph* grVQTimeUnfiltered = DownmixQuadrature(grInputVoltageTemp, lo);

  std::cout<<"First downsampling..."<<std::endl;
  TGraph* grVITimeFirstSample = SampleWaveform(grVITimeUnfiltered, 10*sampleRate);
  TGraph* grVQTimeFirstSample = SampleWaveform(grVQTimeUnfiltered, 10*sampleRate);

  delete grVITimeUnfiltered;
  delete grVQTimeUnfiltered;
  
  // Now we need to filter and then sample these signals
  std::cout<<"Filtering.."<<std::endl;
  TGraph* grVITimeUnsampled = BandPassFilter(grVITimeFirstSample, 0.0, sampleRate/2.0);
  TGraph* grVQTimeUnsampled = BandPassFilter(grVQTimeFirstSample, 0.0, sampleRate/2.0);

  delete grVITimeFirstSample;
  delete grVQTimeFirstSample;
  
  // Now do sampling
  // Use simple linear interpolation for the job
  std::cout<<"Sampling..."<<std::endl;
  grVITime = SampleWaveform(grVITimeUnsampled);
  grVQTime = SampleWaveform(grVQTimeUnsampled);

  delete grVITimeUnsampled;
  delete grVQTimeUnsampled; 
  
  std::cout<<"Adding noise..."<<std::endl;
  AddGaussianNoise(grVITime, noiseTerms);
  AddGaussianNoise(grVQTime, noiseTerms); 
}

rad::Signal::Signal(InducedVoltage iv, LocalOscillator lo, double srate,
		    std::vector<GaussianNoise> noiseTerms, double maxTime) {
  sampleRate = srate;

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
    iv.ResetVoltage();
    iv.GenerateVoltage(lastChunk, thisChunk);
    TGraph* grInputVoltageTemp = iv.GetVoltageGraph();

    std::cout<<"Performing the downmixing..."<<std::endl;
    TGraph* grVITimeUnfiltered = DownmixInPhase(grInputVoltageTemp, lo);
    TGraph* grVQTimeUnfiltered = DownmixQuadrature(grInputVoltageTemp, lo);
  
    std::cout<<"First downsampling"<<std::endl;
    TGraph* grVITimeFirstSample = SampleWaveform(grVITimeUnfiltered, 10*sampleRate, this10Sample);
    TGraph* grVQTimeFirstSample = SampleWaveform(grVQTimeUnfiltered, 10*sampleRate, this10Sample);
    delete grVITimeUnfiltered;
    delete grVQTimeUnfiltered;
    this10Sample = grVITimeFirstSample->GetPointX(grVITimeFirstSample->GetN()-1) + 1/(10*sampleRate);
  
    // Now we need to filter and then sample these signals
    std::cout<<"Filtering.."<<std::endl;
    TGraph* grVITimeUnsampled = BandPassFilter(grVITimeFirstSample, 0.0, sampleRate/2.0);
    TGraph* grVQTimeUnsampled = BandPassFilter(grVQTimeFirstSample, 0.0, sampleRate/2.0);

    delete grVITimeFirstSample;
    delete grVQTimeFirstSample;
  
    // Now do sampling  
    // Use simple linear interpolation for the job
    std::cout<<"Sampling..."<<std::endl;
    TGraph* grVITimeTemp = SampleWaveform(grVITimeUnsampled, sampleRate, thisSample);
    TGraph* grVQTimeTemp = SampleWaveform(grVQTimeUnsampled, sampleRate, thisSample);
    delete grVITimeUnsampled;
    delete grVQTimeUnsampled;
    thisSample = grVITimeTemp->GetPointX(grVITimeTemp->GetN()-1) + 1/sampleRate;

    // Now add the information from these temporary graphs to the larger ones
    for (int i = 0; i < grVITimeTemp->GetN(); i++) {
      grVITime->SetPoint(grVITime->GetN(), grVITimeTemp->GetPointX(i), grVITimeTemp->GetPointY(i));
      grVQTime->SetPoint(grVQTime->GetN(), grVQTimeTemp->GetPointX(i), grVQTimeTemp->GetPointY(i));
    }
    
    delete grVITimeTemp;
    delete grVQTimeTemp;
    
    lastChunk = thisChunk;
    thisChunk += chunkSize;
    if (thisChunk > maxTime) thisChunk = maxTime;
  }
  
  std::cout<<"Adding noise..."<<std::endl;
  AddGaussianNoise(grVITime, noiseTerms);
  AddGaussianNoise(grVQTime, noiseTerms); 
}

rad::Signal::Signal(const Signal &s1) {
  sampleRate = s1.sampleRate;
  grVITime = (TGraph*)s1.grVITime->Clone();
  grVQTime = (TGraph*)s1.grVQTime->Clone();
}

TGraph* rad::Signal::GetVITimeDomain(int firstPoint, int lastPoint) {
  if (firstPoint < 0) firstPoint = 0;
  if (lastPoint < 0 ) lastPoint = grVITime->GetN() - 1;
  
  TGraph* gr = new TGraph();
  setGraphAttr(gr);
  gr->GetXaxis()->SetTitle("Time [s]");
  gr->GetYaxis()->SetTitle("V_{I} [V]");

  for (int i = firstPoint; i <= lastPoint; i++) {
    gr->SetPoint(gr->GetN(), grVITime->GetPointX(i), grVITime->GetPointY(i));
  }
  return gr;
}

TGraph* rad::Signal::GetVQTimeDomain(int firstPoint, int lastPoint) {
  if (firstPoint < 0) firstPoint = 0;
  if (lastPoint < 0 ) lastPoint = grVQTime->GetN() - 1;
  
  TGraph* gr = new TGraph();
  setGraphAttr(gr);
  gr->GetXaxis()->SetTitle("Time [s]");
  gr->GetYaxis()->SetTitle("V_{Q} [V]");

  for (int i = firstPoint; i <= lastPoint; i++) {
    gr->SetPoint(gr->GetN(), grVQTime->GetPointX(i), grVQTime->GetPointY(i));
  }
  return gr;
}

TGraph* rad::Signal::DownmixInPhase(TGraph* grInput, LocalOscillator lo) {
  TGraph* grOut = new TGraph();
  for (int i = 0; i < grInput->GetN(); i++) {
    grOut->SetPoint(i, grInput->GetPointX(i), grInput->GetPointY(i)*lo.GetInPhaseComponent(grInput->GetPointX(i)));
  }
  return grOut;
}

TGraph* rad::Signal::DownmixQuadrature(TGraph* grInput, LocalOscillator lo) {
  TGraph* grOut = new TGraph();
  for (int i = 0; i < grInput->GetN(); i++) {
    grOut->SetPoint(i, grInput->GetPointX(i), grInput->GetPointY(i)*lo.GetQuadratureComponent(grInput->GetPointX(i)));
  }
  return grOut;
}

TGraph* rad::Signal::SampleWaveform(TGraph* grInput) {
  TGraph* grOut = new TGraph();
  double sampleSpacing = 1.0 / sampleRate;
  double sampleTime = grInput->GetPointX(0);

  for (int i = 0; i < grInput->GetN(); i++) {
    double time = grInput->GetPointX(i);
    if (time < sampleTime) continue;
    else if (i == 0) {
      double calcV = grInput->GetPointY(0);
      grOut->SetPoint(grOut->GetN(), sampleTime, calcV);
      sampleTime += sampleSpacing;
    }
    else {
      // Sample the distribution using linear interpolation
      double calcV = grInput->GetPointY(i-1) + (sampleTime - grInput->GetPointX(i-1)) * (grInput->GetPointY(i) - grInput->GetPointY(i-1)) / (time - grInput->GetPointX(i-1));
      grOut->SetPoint(grOut->GetN(), sampleTime, calcV);
      sampleTime += sampleSpacing;      
    }
  }
  
  return grOut;
}

TGraph* rad::Signal::SampleWaveform(TGraph* grInput, const double sRate) {
  TGraph* grOut = new TGraph();
  double sampleSpacing = 1.0 / sRate;
  double sampleTime = grInput->GetPointX(0);
  
  for (int i = 0; i < grInput->GetN(); i++) {
    double time = grInput->GetPointX(i);
    if (time < sampleTime) continue;
    else if (i == 0) {
      double calcV = grInput->GetPointY(0);
      grOut->SetPoint(grOut->GetN(), sampleTime, calcV);
      sampleTime += sampleSpacing;
    }
    else {
      // Sample the distribution using linear interpolation
      double calcV = grInput->GetPointY(i-1) + (sampleTime - grInput->GetPointX(i-1)) * (grInput->GetPointY(i) - grInput->GetPointY(i-1)) / (time - grInput->GetPointX(i-1));
      grOut->SetPoint(grOut->GetN(), sampleTime, calcV);
      sampleTime += sampleSpacing;      
    }
  }

  return grOut;
}

TGraph* rad::Signal::SampleWaveform(TGraph* grInput, const double sRate, const double firstSampleTime) {
  TGraph* grOut = new TGraph();
  double sampleSpacing = 1.0 / sRate;
  double sampleTime = firstSampleTime;
 
  for (int i = 0; i < grInput->GetN(); i++) {
    double time = grInput->GetPointX(i);
    if (time < sampleTime) continue;
    else if (i == 0) {
      double calcV = grInput->GetPointY(0);
      grOut->SetPoint(grOut->GetN(), sampleTime, calcV);
      sampleTime += sampleSpacing;
    }
    else {
      // Sample the distribution using linear interpolation
      double calcV = grInput->GetPointY(i-1) + (sampleTime - grInput->GetPointX(i-1)) * (grInput->GetPointY(i) - grInput->GetPointY(i-1)) / (time - grInput->GetPointX(i-1));
      grOut->SetPoint(grOut->GetN(), sampleTime, calcV);
      sampleTime += sampleSpacing;      
    }
  }

  return grOut;
}

void rad::Signal::AddGaussianNoise(TGraph* grInput, std::vector<GaussianNoise> noiseTerms) {
  for (int i = 0; i < grInput->GetN(); i++) {
    double voltage = grInput->GetPointY(i);
    for (int noise = 0; noise < noiseTerms.size(); noise++) {
      voltage += (noiseTerms.at(noise)).GetNoiseVoltage();  
    }
    grInput->SetPointY(i, voltage);
  }
}

TGraph* rad::Signal::GetVIPowerNorm(const double loadResistance, int firstPoint, int lastPoint) {
  TGraph* grTime = GetVITimeDomain(firstPoint, lastPoint);
  TGraph* grOut = MakePowerSpectrumNorm(grTime);
  delete grTime;
  for (int i = 0; i < grOut->GetN(); i++) {
    grOut->SetPointY(i, grOut->GetPointY(i) / loadResistance);
  }
  setGraphAttr(grOut);
  grOut->GetXaxis()->SetTitle("Frequency [Hz]");
  grOut->GetYaxis()->SetTitle("#frac{V_{I}^{2}}{R} #times (#Deltat)^{2} [W s^{2}]");
  return grOut;
}

TGraph* rad::Signal::GetVQPowerNorm(const double loadResistance, int firstPoint, int lastPoint) {
  TGraph* grTime = GetVQTimeDomain(firstPoint, lastPoint);
  TGraph* grOut = MakePowerSpectrumNorm(grTime);
  delete grTime;
  for (int i = 0; i < grOut->GetN(); i++) {
    grOut->SetPointY(i, grOut->GetPointY(i) / loadResistance);
  }
  setGraphAttr(grOut);
  grOut->GetXaxis()->SetTitle("Frequency [Hz]");
  grOut->GetYaxis()->SetTitle("#frac{V_{Q}^{2}}{R} #times (#Deltat)^{2} [W s^{2}]");
  return grOut;
}

TH2D* rad::Signal::GetVISpectrogram(const double loadResistance, const int NSamplesPerTimeBin) {
  // Get the number of time bins given the specified NSamplesPerTimeBin
  const int nTimeBins = int(floor(double(grVITime->GetN()) / double(NSamplesPerTimeBin)));
  const double firstTime = grVITime->GetPointX(0);
  const double lastTime  = grVITime->GetPointX(nTimeBins * NSamplesPerTimeBin - 1);
  const int nFreqBins = (NSamplesPerTimeBin/2)+1;
  const double deltaF = 1.0 / ((1.0/sampleRate) * NSamplesPerTimeBin);
  const double lastFreq = nFreqBins * deltaF;

  TH2D* h2 = new TH2D("h2", "Spectrogram V_{I}; Time [s]; Frequency [Hz]; #frac{V_{I}^{2}}{R} #times (#Delta t)^{2} [W s^{2}]", nTimeBins, firstTime, lastTime, nFreqBins, 0, lastFreq);
  SetHistAttr(h2);
  // Loop through the time bins and generate the power spectrum at each point
  for (int t = 1; t <= nTimeBins; t++) {
    TGraph* grPower = GetVIPowerNorm(loadResistance, (t-1)*NSamplesPerTimeBin, t*NSamplesPerTimeBin-1);
    // Add this data to the histogram
    for (int i = 0; i < grPower->GetN(); i++) {
      h2->SetBinContent(t, i+1, grPower->GetPointY(i));
    }
    delete grPower;
  }
  
  return h2;
}

TH2D* rad::Signal::GetVQSpectrogram(const double loadResistance, const int NSamplesPerTimeBin) {
  // Get the number of time bins given the specified NSamplesPerTimeBin
  const int nTimeBins = int(floor(double(grVQTime->GetN()) / double(NSamplesPerTimeBin)));
  const double firstTime = grVQTime->GetPointX(0);
  const double lastTime  = grVQTime->GetPointX(nTimeBins * NSamplesPerTimeBin - 1);
  const int nFreqBins = (NSamplesPerTimeBin/2)+1;
  const double deltaF = 1.0 / ((1.0/sampleRate) * NSamplesPerTimeBin);
  const double lastFreq = nFreqBins * deltaF;

  TH2D* h2 = new TH2D("h2", "Spectrogram V_{Q}; Time [s]; Frequency [Hz]; #frac{V_{Q}^{2}}{R} #times (#Delta t)^{2} [W s^{2}]", nTimeBins, firstTime, lastTime, nFreqBins, 0, lastFreq);
  SetHistAttr(h2);
  // Loop through the time bins and generate the power spectrum at each point
  for (int t = 1; t <= nTimeBins; t++) {
    TGraph* grPower = GetVQPowerNorm(loadResistance, (t-1)*NSamplesPerTimeBin, t*NSamplesPerTimeBin-1);
    // Add this data to the histogram
    for (int i = 0; i < grPower->GetN(); i++) {
      h2->SetBinContent(t, i+1, grPower->GetPointY(i));
    }
    delete grPower;
  }
  
  return h2;
}

