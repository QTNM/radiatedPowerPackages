// signalSum.cxx

#include "Antennas/HertzianDipole.h"
#include "BasicFunctions/BasicFunctions.h"
#include "SignalProcessing/LocalOscillator.h"
#include "SignalProcessing/NoiseFunc.h"
#include "SignalProcessing/Signal.h"
#include "TFile.h"
#include "TGraph.h"
#include "TMath.h"
#include "TVector3.h"

using namespace rad;

int main() {
  TVector3 antennaPoint1(0.05, 0.0, 0.0);
  TVector3 dipoleDirZ1(0.0, 1.0, 0.0);
  TVector3 dipoleDirX1(1.0, 0, 0.0);
  HertzianDipole* antenna1 =
      new HertzianDipole(antennaPoint1, dipoleDirX1, dipoleDirZ1, 27.01e9);

  const double antennaAngle2 = 10.0 * TMath::Pi() / 180;
  TVector3 antennaPoint2(0.05 * TMath::Cos(antennaAngle2),
                         0.05 * TMath::Sin(antennaAngle2), 0.0);
  TVector3 dipoleDirZ2(TMath::Sin(antennaAngle2), TMath::Cos(antennaAngle2),
                       0.0);
  TVector3 dipoleDirX2(TMath::Cos(antennaAngle2),
                       -1 * TMath::Sin(antennaAngle2), 0.0);
  HertzianDipole* antenna2 =
      new HertzianDipole(antennaPoint2, dipoleDirX2, dipoleDirZ2, 27.01e9);

  const double antennaAngle3 = -10.0 * TMath::Pi() / 180;
  TVector3 antennaPoint3(0.05 * TMath::Cos(antennaAngle3),
                         0.05 * TMath::Sin(antennaAngle3), 0.0);
  TVector3 dipoleDirZ3(TMath::Sin(antennaAngle3), TMath::Cos(antennaAngle3),
                       0.0);
  TVector3 dipoleDirX3(TMath::Cos(antennaAngle3),
                       -1 * TMath::Sin(antennaAngle3), 0.0);
  HertzianDipole* antenna3 =
      new HertzianDipole(antennaPoint3, dipoleDirX3, dipoleDirZ3, 27.01e9);

  const double loadResistance = 70.0;
  const double noiseTemp = 0.0005;
  LocalOscillator myLO(26.75e9 * 2 * TMath::Pi());
  GaussianNoise noise1(noiseTemp, loadResistance);
  std::vector<GaussianNoise> noiseTerms;
  noiseTerms.push_back(noise1);
  const double sampleRate = 0.75e9;  // Hz

  TString trackFile{"/home/sjones/work/qtnm/trajectories/90DegOnAxis.root"};

  TFile* fout = new TFile("signalSumOutput.root", "recreate");

  Signal mySignal(trackFile, antenna1, myLO, sampleRate, noiseTerms);
  fout->cd();
  TGraph* grVI = mySignal.GetVITimeDomain();
  TGraph* grVQ = mySignal.GetVQTimeDomain();
  TGraph* grVISpec = MakePowerSpectrumNorm(grVI);
  TGraph* grVQSpec = MakePowerSpectrumNorm(grVQ);

  Signal mySignalBoth(trackFile, {antenna1, antenna2}, myLO, sampleRate,
                      noiseTerms);
  fout->cd();
  TGraph* grVIBoth = mySignalBoth.GetVITimeDomain();
  TGraph* grVQBoth = mySignalBoth.GetVQTimeDomain();
  TGraph* grVIBothSpec = MakePowerSpectrumNorm(grVIBoth);
  TGraph* grVQBothSpec = MakePowerSpectrumNorm(grVQBoth);

  Signal mySignalAll(trackFile, {antenna1, antenna2, antenna3}, myLO,
                     sampleRate, noiseTerms);
  fout->cd();
  TGraph* grVIAll = mySignalAll.GetVITimeDomain();
  TGraph* grVQAll = mySignalAll.GetVQTimeDomain();
  TGraph* grVIAllSpec = MakePowerSpectrumNorm(grVIAll);
  TGraph* grVQAllSpec = MakePowerSpectrumNorm(grVQAll);

  grVI->Write("grVI");
  grVQ->Write("grVQ");
  grVISpec->Write("grVISpec");
  grVQSpec->Write("grVQSpec");

  grVIBoth->Write("grVIBoth");
  grVQBoth->Write("grVQBoth");
  grVIBothSpec->Write("grVIBothSpec");
  grVQBothSpec->Write("grVQBothSpec");

  grVIAll->Write("grVIAll");
  grVQAll->Write("grVQAll");
  grVIAllSpec->Write("grVIAllSpec");
  grVQAllSpec->Write("grVQAllSpec");

  fout->Close();
  delete fout;
  return 0;
}
