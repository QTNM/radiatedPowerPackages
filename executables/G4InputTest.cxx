/*
  G4InputTest.cxx

  Check what the results of using G4 input are
*/

#include <getopt.h>

#include <iostream>
#include <string>
#include <vector>

#include "BasicFunctions/BasicFunctions.h"
#include "BasicFunctions/Constants.h"
#include "ElectronDynamics/QTNMFields.h"
#include "ElectronDynamics/TrajectoryGen.h"
#include "SignalProcessing/LocalOscillator.h"
#include "SignalProcessing/Signal.h"
#include "TBranch.h"
#include "TFile.h"
#include "TGraph.h"
#include "TString.h"
#include "TTree.h"
#include "Waveguides/CircularWaveguide.h"

using namespace rad;

int main(int argc, char *argv[]) {
  int opt{};
  std::string outputFileStr{" "};
  while ((opt = getopt(argc, argv, ":o:")) != -1) {
    switch (opt) {
      case 'o':
        outputFileStr = optarg;
        std::cout << "Output file is " << outputFileStr << std::endl;
        break;

      case ':':
        std::cout << "Option needs a value\n";
        break;

      case '?':
        std::cout << "Unknown option: " << optopt << std::endl;
        break;
    }
  }

  TString outputFile{outputFileStr};
  TFile fout(outputFile, "recreate");

  // Define a harmonic trap to match what's in the G4 files
  const double bBkg{1};          // tesla
  const double rCoil{20e-3};     // metres
  const double iCoil{100};       // amps
  const double coilZPos{50e-3};  // metres
  /*
  auto field =
      new BathtubField(rCoil, iCoil, -coilZPos, coilZPos, TVector3(0, 0, bBkg));
  */
  auto field = new HarmonicField(rCoil, iCoil, 1);
  TVector3 centralB{field->evaluate_field_at_point(TVector3(0, 0, 0))};
  const double centralBMag{centralB.Mag()};
  // Calculate central frequency
  const double eKE{18.575e3};  // eV
  const double centralF{CalcCyclotronFreq(eKE, centralBMag)};
  std::cout << "Central magnetic field magnitude = " << centralBMag << " T\n";
  std::cout << "Central frequency = " << centralF / 1e9 << " GHz\n";

  // Define local oscillator
  const double sRate{1e9};                    // Hertz
  const double loFreq{centralF - sRate / 4};  // Hertz
  std::cout << "Local oscillator frequency = " << loFreq / 1e9 << " GHz\n";
  LocalOscillator lo(2 * M_PI * loFreq);

  // Define waveguide
  const double wgLength{12e-2};  // metres
  const double wgRadius{5e-3};   // metres
  auto wg = new CircularWaveguide(wgRadius, wgLength, TVector3(0, 0, 20e-2));

  // Open the G4 file to get the trajectory
  // Define the G4 file path
  TString fileG4{"/home/sjones/work/qtnm/QTNMSim_tests/harmonicOutput.root"};
  TFile fin(fileG4, "read");
  TTree *tr = (TTree *)fin.Get("ntuple/Signal");
  std::vector<double> *vtime = nullptr;
  std::vector<double> *vpx = nullptr;
  std::vector<double> *vpy = nullptr;
  std::vector<double> *vpz = nullptr;
  std::vector<double> *vbx = nullptr;
  std::vector<double> *vby = nullptr;
  std::vector<double> *vbz = nullptr;
  TBranch *bvtime = nullptr;
  TBranch *bvpx = nullptr;
  TBranch *bvpy = nullptr;
  TBranch *bvpz = nullptr;
  TBranch *bvbx = nullptr;
  TBranch *bvby = nullptr;
  TBranch *bvbz = nullptr;
  tr->SetBranchAddress("TimeVec", &vtime, &bvtime);
  tr->SetBranchAddress("PosxVec", &vpx, &bvpx);
  tr->SetBranchAddress("PosyVec", &vpy, &bvpy);
  tr->SetBranchAddress("PoszVec", &vpz, &bvpz);
  tr->SetBranchAddress("BetaxVec", &vbx, &bvbx);
  tr->SetBranchAddress("BetayVec", &vby, &bvby);
  tr->SetBranchAddress("BetazVec", &vbz, &bvbz);
  tr->GetEntry(0);

  TVector3 x0(vpx->at(0), vpy->at(0), vpz->at(0));
  TVector3 v0(vbx->at(0) * TMath::C(), vby->at(0) * TMath::C(),
              vbz->at(0) * TMath::C());
  std::cout << "Start pos = " << x0.X() << ", " << x0.Y() << ", " << x0.Z()
            << std::endl;
  std::cout << "Start vel = " << v0.X() << ", " << v0.Y() << ", " << v0.Z()
            << std::endl;

  auto grX_g4{new TGraph()};
  setGraphAttr(grX_g4);
  grX_g4->SetTitle("G4; Time [s]; X [m]");
  auto grZ_g4{new TGraph()};
  setGraphAttr(grZ_g4);
  grZ_g4->SetTitle("G4; Time [s]; Z [m]");

  const double tAcq{2e-6};  // seconds

  for (size_t g4Ind{0}; g4Ind < vpx->size(); g4Ind++) {
    if (vtime->at(g4Ind) * 1e-9 > tAcq) break;
    grX_g4->SetPoint(grX_g4->GetN(), vtime->at(g4Ind) * 1e-9, vpx->at(g4Ind));
    grZ_g4->SetPoint(grZ_g4->GetN(), vtime->at(g4Ind) * 1e-9, vpz->at(g4Ind));
  }

  fin.Close();

  // Generate a trajectory
  TString trackFile{"/home/sjones/work/qtnm/outputs/G4InputTest/track.root"};
  ElectronTrajectoryGen traj(trackFile, field, x0, v0, 2e-12, tAcq);
  // Now open this file and get the x and y position
  auto grX_h{new TGraph()};
  setGraphAttr(grX_h);
  grX_h->SetTitle("Homemade; Time [s]; X [m]");
  grX_h->SetLineColor(kRed);
  auto grZ_h{new TGraph()};
  setGraphAttr(grZ_h);
  grZ_h->SetTitle("Homemade; Time [s]; Z [m]");
  grZ_h->SetLineColor(kRed);
  TFile fH(trackFile, "read");
  TTree *trH = (TTree*)fH.Get("tree");
  double time{};
  double xPos{}, yPos{}, zPos{};
  trH->SetBranchAddress("time", &time);
  trH->SetBranchAddress("xPos", &xPos);
  trH->SetBranchAddress("yPos", &yPos);
  trH->SetBranchAddress("zPos", &zPos);
  for (int iE{0}; iE < trH->GetEntries(); iE++) {
    trH->GetEntry(iE);
    grX_h->SetPoint(grX_h->GetN(), time, xPos);
    grZ_h->SetPoint(grZ_h->GetN(), time, zPos);
  }
  fH.Close();

  fout.cd();
  grX_g4->Write("grX_g4");
  grX_h->Write("grX_h");
  grZ_g4->Write("grZ_g4");
  grZ_h->Write("grZ_h");

  // Now launch the signal processing
  std::cout << "Starting signal processing with homemade file\n";
  Signal sig1(trackFile, wg, lo, sRate, {}, tAcq);
  std::cout << "Finished signal processing homemade file\n";
  auto grV1{sig1.GetVITimeDomain()};
  auto grVSpec1{sig1.GetVIPowerPeriodogram(1)};
  fout.cd();
  grV1->Write("grV1");
  grVSpec1->Write("grVSpec1");

  std::cout << "Starting signal processing with G4 file\n";
  Signal sig2(fileG4, wg, lo, sRate, {}, tAcq);
  std::cout << "Finished signal processing G4 file\n";

  auto grV2{sig2.GetVITimeDomain()};
  auto grVSpec2{sig2.GetVIPowerPeriodogram(1)};
  fout.cd();
  grV2->Write("grV2");
  grVSpec2->Write("grVSpec2");

  fout.Close();
  return 0;
}