// perpendicularAntenna.cxx
/*
  Tests the possibility of using two perpendicular antennas along with cross-correlation techniques 
  to pick out the signal from the noise
*/

#include <iostream>
#include <vector>

#include "TFile.h"
#include "TAxis.h"
#include "TGraph.h"
#include "TMath.h"
#include "TString.h"
#include "TArrow.h"

#include "FFTtools.h"

#include "SignalProcessing/Signal.h"
#include "SignalProcessing/NoiseFunc.h"
#include "SignalProcessing/LocalOscillator.h"
#include "BasicFunctions/BasicFunctions.h"
#include "Antennas/HalfWaveDipole.h"

using namespace rad;

int main(int argc, char** argv)
{
  std::cout<<"Have "<<argc<<" arguments"<<std::endl;
  for (int i = 0; i < argc; ++i) {
    std::cout<<argv[i]<<std::endl;
  }
  TString outputFile = argv[1];

  const double antennaRadius = 0.02;
  const double antennaAngle1 = 0*TMath::Pi()/180;
  const double antennaAngle2 = 90*TMath::Pi()/180;

  TVector3 antennaPoint1(antennaRadius*TMath::Cos(antennaAngle1), antennaRadius*TMath::Sin(antennaAngle1), 0.0);
  TVector3 antennaDirZ1(-1*TMath::Sin(antennaAngle1), TMath::Cos(antennaAngle1), 0.0);
  TVector3 antennaDirX1(TMath::Cos(antennaAngle1), TMath::Sin(antennaAngle1), 0.0);
  HalfWaveDipole* antenna1 = new HalfWaveDipole(antennaPoint1, antennaDirX1, antennaDirZ1, 27.01e9);
  
  TVector3 antennaPoint2(antennaRadius*TMath::Cos(antennaAngle2), antennaRadius*TMath::Sin(antennaAngle2), 0.0);
  TVector3 antennaDirZ2(-1*TMath::Sin(antennaAngle2), TMath::Cos(antennaAngle2), 0.0);
  TVector3 antennaDirX2(TMath::Cos(antennaAngle2), TMath::Sin(antennaAngle2), 0.0);
  HalfWaveDipole* antenna2 = new HalfWaveDipole(antennaPoint2, antennaDirX2, antennaDirZ2, 27.01e9);

  const double loadResistance = 70.0;
  const double noiseTemp = 4.0;
  LocalOscillator myLO(26.75e9 * 2 * TMath::Pi());
  GaussianNoise noise1(noiseTemp, loadResistance);
  std::vector<GaussianNoise> noiseTerms;
  noiseTerms.push_back(noise1);
  const double sampleRate = 0.75e9; // Hz
  const double tAcq = 54e-6;

  InducedVoltage iv1("/home/sjones/work/qtnm/trajectories/electronTraj60us90Deg.root", antenna1);
  InducedVoltage iv2("/home/sjones/work/qtnm/trajectories/electronTraj60us90Deg.root", antenna2);

  Signal signal1(iv1, myLO, sampleRate, noiseTerms, tAcq);
  Signal signal2(iv2, myLO, sampleRate, noiseTerms, tAcq);

  Signal signal1NoNoise(iv1, myLO, sampleRate, {}, tAcq);
  Signal signal2NoNoise(iv2, myLO, sampleRate, {}, tAcq);
  
  TFile* fout = new TFile(outputFile, "RECREATE");
    
  TGraph* grVI1 = signal1.GetVITimeDomain();
  TGraph* grVI2 = signal2.GetVITimeDomain();
  TGraph* grVI1Spec = signal1.GetVIPowerPeriodogram(loadResistance);
  TGraph* grVI2Spec = signal2.GetVIPowerPeriodogram(loadResistance);

  TGraph* grVI1NoNoise = signal1NoNoise.GetVITimeDomain();
  TGraph* grVI2NoNoise = signal2NoNoise.GetVITimeDomain();
  TGraph* grVI1SpecNoNoise = signal1NoNoise.GetVIPowerPeriodogram(loadResistance);
  TGraph* grVI2SpecNoNoise = signal2NoNoise.GetVIPowerPeriodogram(loadResistance);
  
  fout->cd();
  grVI1->Write("grVI1");
  grVI2->Write("grVI2");
  grVI1Spec->Write("grVI1Spec");
  grVI2Spec->Write("grVI2Spec");

  grVI1NoNoise->Write("grVI1NoNoise");
  grVI2NoNoise->Write("grVI2NoNoise");
  grVI1SpecNoNoise->Write("grVI1SpecNoNoise");
  grVI2SpecNoNoise->Write("grVI2SpecNoNoise");

  TGraph* grCorrelation = FFTtools::getCorrelationGraph(grVI1, grVI2);
  TGraph* grNormCorrelation = FFTtools::getNormalisedCorrelationGraph(grVI1, grVI2);
  TGraph* grNormCorrelationTime = FFTtools::getNormalisedCorrelationGraphTimeDomain(grVI1, grVI2);
  setGraphAttr(grCorrelation);
  setGraphAttr(grNormCorrelation);
  setGraphAttr(grNormCorrelationTime);
  grCorrelation->Write("grCorrelation");
  grNormCorrelation->Write("grNormCorrelation");
  grNormCorrelationTime->Write("grNormCorrelationTime");

  TGraph* grCorrelationFFT = MakeFFTMagGraph(grCorrelation);
  grCorrelationFFT->Write("grCorrelationFFT");
  
  TGraph* grCorrelationNoNoise = FFTtools::getCorrelationGraph(grVI1NoNoise, grVI2NoNoise);
  TGraph* grCorrelationNoNoiseFFT = MakeFFTMagGraph(grCorrelationNoNoise);
  grCorrelationNoNoiseFFT->Write("grCorrelationNoNoiseFFT");
  
  TArrow* ar = new TArrow(27.01e9 - 26.75e9, 19.5e-27, 27.01e9 - 26.75e9, 17e-27, 0.03, "|>");
  ar->SetArrowSize(0.03);
  ar->SetLineWidth(2);
  ar->SetLineColor(kRed);
  ar->SetFillColor(kRed);
    
  fout->Close();
  delete fout;
    
  return 0;
}
