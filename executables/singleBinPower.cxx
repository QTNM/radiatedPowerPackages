// singleBinPower.cxx

#include <vector>

#include "physics/Antennas/HalfWaveDipole.h"
#include "utilities/BasicFunctions/BasicFunctions.h"
#include "physics/SignalProcessing/InducedVoltage.h"
#include "physics/SignalProcessing/LocalOscillator.h"
#include "physics/SignalProcessing/NoiseFunc.h"
#include "physics/SignalProcessing/Signal.h"
#include "TAxis.h"
#include "TFile.h"
#include "TGraph.h"
#include "TString.h"

using namespace rad;

int main(int argc, char** argv) {
  std::cout << "Have " << argc << " arguments" << std::endl;
  for (int i = 0; i < argc; ++i) {
    std::cout << argv[i] << std::endl;
  }
  TString outputFile = argv[1];
  TString trackFile{
      "/home/sjones/work/qtnm/trajectories/electronTraj60us90Deg.root"};

  const double antennaLowerBandwidth{26.95e9};
  const double antennaUpperBandwidth{27.05e9};

  const double acquisitionTime{5e-6};

  TVector3 antennaPoint1(0.02, 0, 0);
  TVector3 antennaDirZ1(0, 1, 0);
  TVector3 antennaDirX1(1, 0, 0);
  HalfWaveDipole* antenna1 =
      new HalfWaveDipole(antennaPoint1, antennaDirX1, antennaDirZ1, 27.01e9);
  antenna1->SetBandwidth(antennaLowerBandwidth, antennaUpperBandwidth);

  LocalOscillator myLO(26.75e9 * 2 * PI);
  const double loadResistance = 70.0;
  const double sampleRate = 750e6;  // Hz
  const double noiseTemp = 4.0;
  GaussianNoise noise1(noiseTemp, loadResistance);
  std::vector<GaussianNoise> noiseTerms;
  noiseTerms.push_back(noise1);

  Signal signal(trackFile, antenna1, myLO, sampleRate, noiseTerms);
  Signal signalNoNoise(trackFile, antenna1, myLO, sampleRate, {});

  TFile* fout = new TFile(outputFile, "RECREATE");
  fout->cd();

  TGraph* grVIPeriodogram = signal.GetVIPowerPeriodogram(loadResistance);
  TGraph* grVIPeriodogramNoNoise =
      signalNoNoise.GetVIPowerPeriodogram(loadResistance);
  std::cout << "Power integral " << SumPower(grVIPeriodogram) << std::endl;
  std::cout << "Power integral no noise " << SumPower(grVIPeriodogramNoNoise)
            << std::endl;
  grVIPeriodogram->Write("grVIPeriodogram");
  grVIPeriodogramNoNoise->Write("grVIPeriodogramNoNoise");

  fout->Close();
  delete fout;

  return 0;
}
