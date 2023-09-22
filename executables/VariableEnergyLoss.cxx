/*
  VariableEnergyLoss.cxx


*/

#include <getopt.h>
#include <unistd.h>

#include <iostream>
#include <string>

#include "BasicFunctions/BasicFunctions.h"
#include "BasicFunctions/Constants.h"
#include "ElectronDynamics/BorisSolver.h"
#include "ElectronDynamics/QTNMFields.h"
#include "ElectronDynamics/TrajectoryGen.h"
#include "TFile.h"
#include "TMath.h"
#include "TString.h"
#include "TTree.h"
#include "TTreeReader.h"
#include "TTreeReaderValue.h"
#include "Waveguides/CircularCavity.h"

using namespace rad;

int main(int argc, char *argv[]) {
  // Input parsing
  int opt{};
  std::string outputStr{" "};
  while ((opt = getopt(argc, argv, ":o:")) != -1) {
    switch (opt) {
      case 'o':
        outputStr = optarg;
        std::cout << "Output file is " << outputStr << std::endl;
        break;

      case ':':
        std::cout << "Option needs a value" << std::endl;
        break;

      case '?':
        std::cout << "Unknown option: " << optopt << std::endl;
        break;
    }
  }

  // Create output file
  TString outputFile{outputStr};
  TFile fout(outputFile, "RECREATE");

  // Generate the field
  // This is chosen so we get the correct frequency
  const double rCoil{3e-2};      // metres
  const double trapDepth{3e-3};  // T
  const double fieldMag{0.7};    // T
  const double iCoil{2 * trapDepth * rCoil / MU0};
  auto field = new HarmonicField(rCoil, iCoil, fieldMag);
  // Get the field at the centre of the trap
  TVector3 centreField{field->evaluate_field_at_point(TVector3(0, 0, 0))};

  // Electron setup
  const double ke{18.6e3};  // eV
  const double speed{GetSpeedFromKE(ke, ME)};
  TVector3 initVel(speed, 0, 0);
  const double gyroradius{GetGyroradius(initVel, centreField, ME)};
  TVector3 initPos(0, gyroradius, 0);

  // Calculate cyclotron frequency
  const double cycFreq{CalcCyclotronFreq(ke, centreField.Mag())};  // Hertz
  std::cout << "Cyclotron frequency = " << cycFreq / 1e9 << " GHz\n";

  // Define the cavity
  const double cavityRadius{5e-3};  // metres
  const double p11Prime{GetBesselPrimeZero(1, 1)};
  const double cavityLength{TMath::Pi() /
                            sqrt(pow(TMath::TwoPi() * cycFreq / TMath::C(), 2) -
                                 pow(p11Prime / cavityRadius, 2))};
  std::cout << "Cavity length = " << cavityLength * 1e3 << " mm\n";
  TVector3 probePosition(0.5 * cavityRadius, 0, 0);
  auto cavity = new CircularCavity(cavityRadius, cavityLength, probePosition);

  const double simulationTime{1e-7};  // seconds
  const double simStepSize{2e-12};    // seconds
  TString trackFileFreeSpace{
      "~/work/qtnm/Cavities/VariableEnergyLoss/trackFreeSpace.root"};
  // Generate free space file
  const clock_t beginFreeSpace{clock()};
  std::cout << "Starting free space generation\n";
  ElectronTrajectoryGen trajFreeSpace(trackFileFreeSpace, field, initPos,
                                      initVel, simStepSize, simulationTime,
                                      true);
  const clock_t endFreeSpace{clock()};
  const double freeSpaceTime{double(endFreeSpace - beginFreeSpace) /
                             CLOCKS_PER_SEC};
  std::cout << "Finished free space generation in " << freeSpaceTime
            << " s\n\n";

  const clock_t beginCavity{clock()};
  TString trackFileCavity{
      "~/work/qtnm/Cavities/VariableEnergyLoss/trackCavity.root"};
  std::cout << "Starting cavity generation\n";
  ElectronTrajectoryGen trajCavity(trackFileCavity, field, initPos, initVel,
                                   simStepSize, simulationTime, true, cavity);
  const clock_t endCavity{clock()};
  const double cavityTime{double(endCavity - beginCavity) / CLOCKS_PER_SEC};
  std::cout << "Finished cavity generation in " << cavityTime << " s\n";

  // Now read in the files and figure out the rate of enery loss
  TFile fFreeSpace(trackFileFreeSpace, "read");
  TTreeReader readFreeSpace("tree", &fFreeSpace);
  TTreeReaderValue<double> tFree(readFreeSpace, "time");
  TTreeReaderValue<double> xFree(readFreeSpace, "xVel");
  TTreeReaderValue<double> yFree(readFreeSpace, "yVel");
  TTreeReaderValue<double> zFree(readFreeSpace, "zVel");

  double E0{ke * TMath::Qe()};

  auto grPowerFree = new TGraph();
  setGraphAttr(grPowerFree);
  grPowerFree->SetTitle("Free space; Time [s]; P_{rad} [fW]");
  while (readFreeSpace.Next()) {
    TVector3 vel(*xFree, *yFree, *zFree);
    double gamma1{1 / sqrt(1 - pow(vel.Mag() / TMath::C(), 2))};
    double E1{(gamma1 - 1) * ME * TMath::C() * TMath::C()};
    double power{(E0 - E1) * 1e15 / simStepSize};
    if (power > 0.01) {
      grPowerFree->SetPoint(grPowerFree->GetN(), *tFree, power);
    }
    E0 = E1;
  }
  fFreeSpace.Close();
  fout.cd();
  grPowerFree->Write("grPowerFree");

  TFile fCavity(trackFileCavity, "read");
  TTreeReader readCavity("tree", &fCavity);
  TTreeReaderValue<double> tCav(readCavity, "time");
  TTreeReaderValue<double> xCav(readCavity, "xVel");
  TTreeReaderValue<double> yCav(readCavity, "yVel");
  TTreeReaderValue<double> zCav(readCavity, "zVel");

  auto grPowerCav = new TGraph();
  setGraphAttr(grPowerCav);
  grPowerCav->SetTitle("Cavity; Time [s]; P_{rad} [fW]");
  grPowerCav->SetLineColor(kRed);
  E0 = ke * TMath::Qe();
  while (readCavity.Next()) {
    TVector3 vel(*xCav, *yCav, *zCav);
    double gamma1{1 / sqrt(1 - pow(vel.Mag() / TMath::C(), 2))};
    double E1{(gamma1 - 1) * ME * TMath::C() * TMath::C()};
    double power{(E0 - E1) * 1e15 / simStepSize};
    if (power > 0.01) {
      grPowerCav->SetPoint(grPowerCav->GetN(), *tCav, power);
    }
    E0 = E1;
  }
  fCavity.Close();
  fout.cd();
  grPowerCav->Write("grPowerCav");

  fout.Close();
  return 0;
}