/*
  NicePenningTraj.cxx

  Make a nice trajectory plot of an electron in a Penning trap
*/

#include <iostream>

#include "BasicFunctions/BasicFunctions.h"
#include "BasicFunctions/Constants.h"
#include "ElectronDynamics/PenningTraps.h"
#include "ElectronDynamics/QTNMFields.h"
#include "ElectronDynamics/TrajectoryGen.h"
#include "TFile.h"
#include "TGraph2D.h"
#include "TMath.h"
#include "TString.h"
#include "TTree.h"
#include "TTreeReader.h"
#include "TTreeReaderValue.h"

using namespace rad;

int main() {
  TString outputPath{"/home/sjones/work/qtnm/outputs/NicePenningTraj/"};
  TFile fout(outputPath + "/plots.root", "recreate");

  // Electron kinematics
  const double eKE{100};                         // eV
  const double eSpeed{GetSpeedFromKE(eKE, ME)};  // m s^-1
  TVector3 initPos(0.02, 0, 0);

  // Set up Penning trap
  const double BFieldMag{0.01};  // T
  const double v0{-100};         // V
  const double z0{0.1};          // m
  const double rho0{0.05};       // m
  auto trap = new IdealPenningTrap(BFieldMag, v0, rho0, z0);

  // Calculate minimum pitch angle
  const double thetaBot{acos(sqrt(-TMath::Qe() * v0 / (eKE * TMath::Qe()))) *
                        180 / TMath::Pi()};  // degrees
  std::cout << "Minimum trap angle = " << thetaBot << " degrees\n";

  const double angle{88.0 * TMath::Pi() / 180};
  TVector3 initVel(eSpeed * sin(angle), 0, eSpeed * cos(angle));
  const double stepSize{1e-12};
  const double simTime{1e-6};
  TString trackFile{outputPath + "/track.root"};
  ElectronTrajectoryGen traj(trackFile, trap, initPos, initVel, stepSize,
                             simTime, false);

  // Now read in track file
  TFile fin(trackFile, "READ");
  // Get the TTree
  TTreeReader reader("tree", &fin);
  TTreeReaderValue<double> time(reader, "time");
  TTreeReaderValue<double> xPos(reader, "xPos");
  TTreeReaderValue<double> yPos(reader, "yPos");
  TTreeReaderValue<double> zPos(reader, "zPos");

  // Use this to draw a nice 3D electron trajectory
  TGraph2D helix;
  helix.SetTitle("; x [cm]; y [cm]; z [cm]");

  const double lastTime{8e-8};  // We don't want to draw too much
  while (reader.Next()) {
    if (*time > lastTime) break;

    helix.SetPoint(helix.GetN(), *xPos * 100, *yPos * 100, *zPos * 100);
  }
  fin.Close();

  fout.cd();
  helix.Write("helix");

  fout.Close();
  return 0;
}