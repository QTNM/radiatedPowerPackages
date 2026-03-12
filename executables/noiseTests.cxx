// noiseTests.cxx

#include "TFile.h"
#include "TGraph.h"
#include "TString.h"
#include "physics/Antennas/HertzianDipole.h"
#include "physics/SignalProcessing/LocalOscillator.h"
#include "physics/SignalProcessing/NoiseFunc.h"
#include "physics/SignalProcessing/Signal.h"
#include "utilities/BasicCore/Constants.h"
#include "utilities/ROOTUtils/FFTAnalysis.h"

using namespace rad;

int main() {
  TVector3 antennaPoint(0.02, 0.0, 0.0);
  TVector3 dipoleDirZ(0.0, 1.0, 0.0);
  TVector3 dipoleDirX(1.0, 0.0, 0.0);
  HertzianDipole* myAntenna =
      new HertzianDipole(antennaPoint, dipoleDirX, dipoleDirZ, 27.01e9);

  const double loadResistance = 70.0;
  const double noiseTemp = 0.001;
  LocalOscillator myLO(26.75e9 * 2 * PI);
  GaussianNoise noise1(noiseTemp, loadResistance);
  std::vector<GaussianNoise> noiseTerms;
  noiseTerms.push_back(noise1);
  const double sampleRate = 0.75e9;  // Hz

  TString trackFile{"/home/sjones/work/qtnm/trajectories/90DegOnAxis.root"};

  TFile* fout = new TFile("noiseTestOutput.root", "RECREATE");
  Signal mySignal(trackFile, myAntenna, myLO, sampleRate, noiseTerms);

  fout->cd();
  TGraph* grVI = mySignal.GetVITimeDomain();
  TGraph* grVQ = mySignal.GetVQTimeDomain();
  TGraph* grVISpec = MakePowerSpectrumNorm(grVI);
  TGraph* grVQSpec = MakePowerSpectrumNorm(grVQ);

  grVI->Write("grVI");
  grVQ->Write("grVQ");
  grVISpec->Write("grVISpec");
  grVQSpec->Write("grVQSpec");

  fout->Close();
  delete fout;

  return 0;
}
