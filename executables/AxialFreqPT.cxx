/*
  AxialFreqPT.cxx

  Check the variation in axial frequency as a function of pitch angle
*/

#include <iostream>
#include <memory>

#include "BasicFunctions/BasicFunctions.h"
#include "BasicFunctions/Constants.h"
#include "ElectronDynamics/PenningTraps.h"
#include "ElectronDynamics/TrajectoryGen.h"
#include "TFile.h"
#include "TGraph.h"
#include "TString.h"
#include "TTree.h"

using namespace rad;

double GetAxialFreq(TString inputFile) {
  auto fin = std::make_unique<TFile>(inputFile, "read");
  auto tr = (TTree *)fin->Get("tree");
  double time{0};
  double zPos{0};
  tr->SetBranchAddress("time", &time);
  tr->SetBranchAddress("zPos", &zPos);

  auto gr = new TGraph();
  for (int n{0}; n < tr->GetEntries(); n++) {
    tr->GetEntry(n);
    gr->SetPoint(n, time, zPos);
  }
  // Fourier transform
  auto grP = MakePowerSpectrumPeriodogram(gr);
  delete gr;

  // Find biggest frequency component
  double maxPower{0};
  double maxPowerFreq{0};
  for (int i{0}; i < grP->GetN(); i++) {
    if (grP->GetPointY(i) > maxPower) {
      maxPower = grP->GetPointY(i);
      maxPowerFreq = grP->GetPointX(i);
    }
  }
  delete grP;

  delete tr;
  fin->Close();
  return maxPowerFreq;
}

TGraph *GetZPosGraph(TString inputFile) {
  auto fin = std::make_unique<TFile>(inputFile, "read");
  auto tr = (TTree *)fin->Get("tree");
  double time{0};
  double zPos{0};
  tr->SetBranchAddress("time", &time);
  tr->SetBranchAddress("zPos", &zPos);

  auto gr = new TGraph();
  setGraphAttr(gr);
  for (int n{0}; n < tr->GetEntries(); n++) {
    tr->GetEntry(n);
    gr->SetPoint(n, time, zPos);
  }
  delete tr;
  fin->Close();
  return gr;
};

TGraph *GetRhoPowerSpec(TString inputFile) {
  auto fin = std::make_unique<TFile>(inputFile, "read");
  auto tr = (TTree *)fin->Get("tree");
  double time{0};
  double xPos{0};
  double yPos{0};
  tr->SetBranchAddress("time", &time);
  tr->SetBranchAddress("xPos", &xPos);
  tr->SetBranchAddress("yPos", &yPos);

  auto gr = new TGraph();
  for (int n{0}; n < tr->GetEntries(); n++) {
    tr->GetEntry(n);
    gr->SetPoint(n, time, sqrt(xPos * xPos + yPos * yPos));
  }
  delete tr;
  fin->Close();

  auto grP = MakePowerSpectrumPeriodogram(gr);
  delete gr;
  return grP;
};

int main(int argc, char *argv[]) {
  TString outputFile{argv[1]};
  auto fout = std::make_unique<TFile>(outputFile, "recreate");

  // Electron kinematics
  const double eKE{18.6e3};                      // eV
  const double eSpeed{GetSpeedFromKE(eKE, ME)};  // m s^-1
  TVector3 initPos(0, 0, 0);

  // Set up Penning trap
  const double BFieldMag{1};  // T
  const double v0{-100};      // V
  const double z0{0.1};       // m
  const double rho0{0.05};    // m
  auto trap = new IdealPenningTrap(BFieldMag, v0, rho0, z0);

  // Calculate minimum pitch angle
  const double thetaBot{acos(sqrt(-TMath::Qe() * v0 / (eKE * TMath::Qe()))) *
                        180 / TMath::Pi()};  // degrees
  std::cout << "Minimum trap angle = " << thetaBot << " degrees\n";

  // Set up angle scan
  const double scanAngleMin{86};
  const double scanAngleMax{90};
  const int nScanPnts{20};
  auto grAngleFreq = new TGraph();
  setGraphAttr(grAngleFreq);
  grAngleFreq->SetTitle(
      "Axial frequency in a Penning trap; #theta [degrees]; f_{z} [MHz]");

  // Scan through angles
  for (int iP{0}; iP < nScanPnts; iP++) {
    double scanAngleDeg{scanAngleMin + (scanAngleMax - scanAngleMin) *
                                           double(iP) / double(nScanPnts - 1)};
    std::cout << "Pitch angle = " << scanAngleDeg << " degrees\n";

    double scanAngle{scanAngleDeg * TMath::Pi() / 180};
    TVector3 initVel(eSpeed * sin(scanAngle), 0, eSpeed * cos(scanAngle));
    const double stepSize{1e-12};
    const double simTime{1e-6};
    TString trackFile{
        "/home/sjones/work/qtnm/outputs/PenningTrapTraj/track.root"};
    ElectronTrajectoryGen traj(trackFile, trap, initPos, initVel, stepSize,
                               simTime);
    double fz{GetAxialFreq(trackFile)};
    grAngleFreq->SetPoint(iP, scanAngleDeg, fz / 1e6);

    auto grZ = GetZPosGraph(trackFile);
    grZ->SetTitle(Form("%.2f degrees; Time [s]; z [m]", scanAngleDeg));

    auto grRhoP = GetRhoPowerSpec(trackFile);
    grRhoP->SetTitle(
        Form("%.2f degrees; Frequency [Hz]; [A. U.]", scanAngleDeg));

    fout->cd();
    grZ->Write(Form("grZ%d", iP));
    grRhoP->Write(Form("grRhoP%d", iP));
    delete grZ;
    delete grRhoP;
  }

  fout->cd();
  grAngleFreq->Write("grAngleFreq");

  fout->Close();
  return 0;
}