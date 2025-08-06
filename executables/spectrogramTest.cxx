// spectrogramTest.cxx

#include "TFile.h"
#include "TGraph.h"

#include "physics/FieldClasses/FieldClasses.h"
#include "physics/SignalProcessing/Spectrogram.h"
#include "physics/SignalProcessing/NoiseFunc.h"
#include "physics/Antennas/HertzianDipole.h"

using namespace rad;

int main() {
  TVector3 antennaPoint(0.02, 0.0, 0.0);
  TVector3 dipoleDirZ(0.0, 1.0, 0.0);
  TVector3 dipoleDirX(1.0, 0.0, 0.0);
  HertzianDipole* myAntenna = new HertzianDipole(antennaPoint, dipoleDirX, dipoleDirZ, 27.01e9);
  
  GaussianNoise* noise1 = new GaussianNoise(0.01, 70.0);
  
  FieldPoint *fp = new FieldPoint("/home/sjones/work/qtnm/trajectories/electronTraj600us_89deg.root", myAntenna);
  fp->GenerateFields(0, 50e-6);
  std::cout<<"Generated the fields for the field point"<<std::endl;
  
  TFile* fout = new TFile("spectrogramOutput.root", "RECREATE");
  std::cout<<"Opened the output file"<<std::endl;

  std::vector<GaussianNoise*> noiseTerms;
  noiseTerms.push_back(noise1);
  Spectrogram spec(fp);
  std::cout<<"Created spectrogram class"<<std::endl;
  TH2D* h2spec = spec.MakeSpectrogram(259072, false);
  std::cout<<"Created spectrogram"<<std::endl;

  fout->cd();
  h2spec->Write("h2spec");
  
  fout->Close();
  delete fout;
  
  return 0;
}
