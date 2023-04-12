/*
  PenningTrapTraj.cxx
*/

#include <memory>

#include "BasicFunctions/BasicFunctions.h"
#include "ElectronDynamics/PenningTraps.h"
#include "ElectronDynamics/TrajectoryGen.h"
#include "TFile.h"
#include "TGraph.h"
#include "TGraph2D.h"
#include "TMath.h"
#include "TString.h"
#include "TTree.h"

using namespace rad;

int main(int argc, char *argv[]) {
  TString outputFile{argv[1]};
  auto fout = std::make_unique<TFile>(outputFile, "recreate");

  // Electron kinematics
  const double eKE{18.6e3};
  const double eSpeed{GetSpeedFromKE(eKE, ME)};
  const double pitchAngleDeg{88};
  const double pitchAngle{pitchAngleDeg * TMath::Pi() / 180};
  TVector3 initVel(eSpeed * sin(pitchAngle), 0, eSpeed * cos(pitchAngle));
  TVector3 initPos(0, 0, 0);

  // Set up Penning trap
  const double BFieldMag{1};  // T
  const double v0{-1e3};      // V
  const double z0{0.1};       // m
  const double rho0{0.05};    // m
  auto trap = new IdealPenningTrap(BFieldMag, v0, rho0, z0);

  // Calculate trap frequencies
  const double omegaZ{sqrt(-TMath::Qe() * v0 / (ME * z0 * z0))};
  const double fZ{omegaZ / (2 * TMath::Pi())};
  const double fCFree{CalcCyclotronFreq(eKE, BFieldMag)};
  const double omegaC{2 * TMath::Pi() * fCFree};
  const double omegaPlus{
      0.5 * (omegaC + sqrt(omegaC * omegaC - 2 * omegaZ * omegaZ))};
  const double fPlus{omegaPlus / (2 * TMath::Pi())};
  const double omegaMinus{
      0.5 * (omegaC - sqrt(omegaC * omegaC - 2 * omegaZ * omegaZ))};
  const double fMinus{omegaMinus / (2 * TMath::Pi())};

  std::cout << "Free cyclotron frequency = " << fCFree / 1e9 << " GHz\n";
  std::cout << "Axial frequency = " << fZ / 1e6 << " MHz\n";
  std::cout << "Modified cyclotron freq = " << fPlus / 1e9 << " GHz\n";
  std::cout << "Magnetron freq = " << fPlus / 1e9 << " GHz\n";

  // Simulation details
  const double stepSize{1e-12};
  const double simTime{5e-8};
  TString trackFile{
      "/home/sjones/work/qtnm/outputs/PenningTrapTraj/track.root"};
  ElectronTrajectoryGen traj(trackFile, trap, initPos, initVel, stepSize,
                             simTime);

  // Plot the z position and the radial position as a function of time
  auto tFile = std::make_unique<TFile>(trackFile, "READ");
  auto tr = (TTree *)tFile->Get("tree");
  double time;
  double xPos, yPos, zPos;
  double xVel, yVel, zVel;
  tr->SetBranchAddress("time", &time);
  tr->SetBranchAddress("xPos", &xPos);
  tr->SetBranchAddress("yPos", &yPos);
  tr->SetBranchAddress("zPos", &zPos);
  tr->SetBranchAddress("xVel", &xVel);
  tr->SetBranchAddress("yVel", &yVel);
  tr->SetBranchAddress("zVel", &zVel);

  auto grZ = new TGraph();
  setGraphAttr(grZ);
  grZ->SetTitle("Z; Time [s]; z [m]");
  auto grRho = new TGraph();
  setGraphAttr(grRho);
  grRho->SetTitle("#rho; Time [s]; #rho [m]");
  auto grHelix = new TGraph2D();
  grHelix->SetTitle("Penning trap trajectory; x [m]; y [m]; z[m]");
  auto grPitch = new TGraph();
  setGraphAttr(grPitch);
  grPitch->SetTitle("Pitch angle; Time [s]; #theta [deg]");

  // KE
  auto grKE = new TGraph();
  setGraphAttr(grKE);
  grKE->SetTitle("KE; Time [s]; KE [eV]");
  auto grKEPar = new TGraph();
  setGraphAttr(grKEPar);
  grKEPar->SetTitle("KE_{par}; Time [s]; KE_{par} [eV]");
  auto grKEPerp = new TGraph();
  setGraphAttr(grKEPerp);
  grKEPerp->SetTitle("KE_{perp}; Time [s]; KE_{perp} [eV]");

  for (unsigned int e{0}; e < tr->GetEntries(); e++) {
    tr->GetEntry(e);
    grZ->SetPoint(e, time, zPos);
    grRho->SetPoint(e, time, sqrt(xPos * xPos + yPos * yPos));
    grPitch->SetPoint(
        e, time,
        atan2(sqrt(xVel * xVel + yVel * yVel), zVel) * 180 / TMath::Pi());
    grHelix->SetPoint(e, xPos, yPos, zPos);

    double speed{sqrt(xVel * xVel + yVel * yVel + zVel * zVel)};
    double gamma{1 / sqrt(1 - pow(speed / TMath::C(), 2))};
    double KE{(gamma - 1) * ME * TMath::C() * TMath::C()};
    double speedPerp{sqrt(xVel * xVel + yVel * yVel)};
    grKE->SetPoint(e, time, KE / TMath::Qe());
    grKEPar->SetPoint(e, time, 0.5 * ME * zVel * zVel / TMath::Qe());
    grKEPerp->SetPoint(e, time, 0.5 * ME * speedPerp * speedPerp / TMath::Qe());
  }
  delete tr;
  tFile->Close();

  fout->cd();
  grZ->Write("grZ");
  grRho->Write("grRho");
  grPitch->Write("grPitch");
  grHelix->Write("grHelix");
  grKE->Write("grKE");
  grKEPar->Write("grKEPar");
  grKEPerp->Write("grKEPerp");

  fout->Close();
  return 0;
}