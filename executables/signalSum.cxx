// signalSum.cxx

#include <cmath>

#include "TFile.h"
#include "TGraph.h"
#include "TVector3.h"
#include "physics/Antennas/HertzianDipole.h"
#include "physics/SignalProcessing/LocalOscillator.h"
#include "physics/SignalProcessing/NoiseFunc.h"
#include "physics/SignalProcessing/Signal.h"
#include "utilities/ROOTUtils/FFTAnalysis.h"

using namespace rad;

int main() {
  TVector3 antennaPoint1(0.05, 0.0, 0.0);
  TVector3 dipoleDirZ1(0.0, 1.0, 0.0);
  TVector3 dipoleDirX1(1.0, 0, 0.0);
  HertzianDipole* antenna1 =
      new HertzianDipole(antennaPoint1, dipoleDirX1, dipoleDirZ1, 27.01e9);

  const double antennaAngle2 = 10.0 * PI / 180;
  TVector3 antennaPoint2(0.05 * std::cos(antennaAngle2),
                         0.05 * std::sin(antennaAngle2), 0.0);
  TVector3 dipoleDirZ2(std::sin(antennaAngle2), std::cos(antennaAngle2), 0.0);
  TVector3 dipoleDirX2(std::cos(antennaAngle2), -1 * std::sin(antennaAngle2),
                       0.0);
  HertzianDipole* antenna2 =
      new HertzianDipole(antennaPoint2, dipoleDirX2, dipoleDirZ2, 27.01e9);

  const double antennaAngle3 = -10.0 * PI / 180;
  TVector3 antennaPoint3(0.05 * std::cos(antennaAngle3),
                         0.05 * std::sin(antennaAngle3), 0.0);
  TVector3 dipoleDirZ3(std::sin(antennaAngle3), std::cos(antennaAngle3), 0.0);
  TVector3 dipoleDirX3(std::cos(antennaAngle3), -1 * std::sin(antennaAngle3),
                       0.0);
  HertzianDipole* antenna3 =
      new HertzianDipole(antennaPoint3, dipoleDirX3, dipoleDirZ3, 27.01e9);

  const double loadResistance = 70.0;
  const double noiseTemp = 0.0005;
  LocalOscillator myLO(26.75e9 * 2 * PI);
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
