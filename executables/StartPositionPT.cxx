/*
  StartPositionPT.cxx

  Look at how changing the start position of our electrons affects the radiated
  frequencies. Are we able to recover this information
*/

#include <memory>

#include "BasicFunctions/BasicFunctions.h"
#include "ElectronDynamics/PenningTraps.h"
#include "ElectronDynamics/TrajectoryGen.h"
#include "TF1.h"
#include "TFile.h"
#include "TGraph.h"
#include "TString.h"
#include "TTree.h"

using namespace rad;

int main() {
  TString outputDir{"~/work/qtnm/PenningTrap/StartPositionPT"};
  TString outputFile{outputDir + "/plots.root"};
  auto fout = std::make_unique<TFile>(outputFile, "recreate");

  // Define the Penning trap
  // Make it really deep so we can see the effects on the energy
  const double BFieldMag{1};  // T
  const double v0{-1e3};      // V
  const double z0{0.15};      // m
  const double rho0{0.05};    // m
  auto trap = new IdealPenningTrap(BFieldMag, v0, rho0, z0);
  const double fZ{(1 / (2 * TMath::Pi())) *
                  sqrt(-TMath::Qe() * v0 / (ME * z0 * z0))};
  std::cout << "Axial frequency = " << fZ / 1e6 << " MHz\n";

  // Electron kinematics
  const double eKE{18.6e3};                      // eV
  const double eSpeed{GetSpeedFromKE(eKE, ME)};  // m s^-1
  const double pitchAngleDeg{89};                // degrees
  const double pitchAngle{pitchAngleDeg * TMath::Pi() / 180};
  // Initial velocity vector
  TVector3 initVel(eSpeed * sin(pitchAngle), 0, -eSpeed * cos(pitchAngle));
  // Gyroradius
  const double rg{GetGyroradius(
      initVel, trap->evaluate_field_at_point(TVector3(0, 0, 0)), ME)};

  // Now scan through various positions and measure the axial frequency
  const double zMax{0.12};
  const double zMin{-0.12};
  const unsigned int nZPnts{2};

  auto grAxFreq = new TGraph();
  setGraphAttr(grAxFreq);
  grAxFreq->SetTitle("Axial frequency variation; z_{start} [cm]; f_{z} [MHz]");
  grAxFreq->SetMarkerStyle(20);

  /*
  for (int iZ{0}; iZ < nZPnts; iZ++) {
    double z{zMin + (zMax - zMin) * double(iZ) / double(nZPnts - 1)};
    TVector3 initPos(0, -rg, z);
    std::cout << "z = " << z * 100 << " cm\n";

    // Generate the trajectory
    const double stepSize{1e-12};
    const double simTime{1e-6};
    TString trackFile{outputDir + Form("/track%d.root", iZ)};
    ElectronTrajectoryGen traj(trackFile, trap, initPos, initVel, stepSize,
                               simTime);

    // Now get the radial position as a function of time
    auto fin = std::make_unique<TFile>(trackFile, "read");
    TTree *tr = (TTree *)fin->Get("tree");
    double time{0};
    double xPos{0}, yPos{0}, zPos{0};
    tr->SetBranchAddress("time", &time);
    tr->SetBranchAddress("xPos", &xPos);
    tr->SetBranchAddress("yPos", &yPos);
    tr->SetBranchAddress("zPos", &zPos);

    auto grZ = new TGraph();
    setGraphAttr(grZ);
    grZ->SetTitle(Form("z_{start} = %.0f cm; Time [s]; z [m]", z * 100));
    auto grRho = new TGraph();
    setGraphAttr(grRho);
    grRho->SetTitle(Form("z_{start} = %.0f cm; Time [s]; #rho [m]", z * 100));
    auto grX = new TGraph();
    setGraphAttr(grX);
    grX->SetTitle(Form("z_{start} = %.0f cm; Time [s]; x [m]", z * 100));

    // Loop through tree entries
    double maxAmp{0};
    for (int iE{0}; iE < tr->GetEntries(); iE++) {
      tr->GetEntry(iE);
      grX->SetPoint(iE, time, xPos);
      grZ->SetPoint(iE, time, zPos);

      if (zPos > maxAmp) maxAmp = zPos;

      double rho{sqrt(xPos * xPos + yPos * yPos)};
      grRho->SetPoint(iE, time, rho);
    }

    auto fAx = new TF1("fAx", "[0] * sin(2 * 3.14159 * [1] * x + [2])");
    fAx->SetParameter(0, maxAmp);
    fAx->SetParameter(1, fZ);
    fAx->SetParameter(2, 1);
    grZ->Fit(fAx);
    grAxFreq->SetPoint(iZ, z * 100, fAx->GetParameter(1));

    fout->cd();
    grX->Write(Form("grX%d", iZ));
    grZ->Write(Form("grZ%d", iZ));
    grRho->Write(Form("grRho%d", iZ));
    delete fAx;
    delete grZ;
    delete grRho;

    auto grXPower{MakePowerSpectrumPeriodogram(grX)};
    delete grX;
    grXPower->SetTitle(
        Form("z_{start} = %.0f cm; f [Hz]; FFT magnitude [A. U.]", z * 100));

    fout->cd();
    grXPower->Write(Form("grXPower%d", iZ));
    delete grXPower;

    delete tr;
    fin->Close();
  }
  fout->cd();
  grAxFreq->Write("grAxFreq");
  */
  //////////////////////////////////////////////////////////////////////////////
  /// Generate two more trajectories to show how two electrons with different
  /// initial energies can end up with the same trajectory
  //////////////////////////////////////////////////////////////////////////////
  const double theta{88 * TMath::Pi() / 180};
  TVector3 startVel88(eSpeed * sin(theta), 0, -eSpeed * cos(theta) * 1.0339054);
  double ke88{(1 / sqrt(1 - pow(eSpeed / TMath::C(), 2)) - 1) * ME *
              TMath::C() * TMath::C()};
  TVector3 startVel90(eSpeed * sin(theta) * 1.000041899, 0, 0);
  double ke90{(1 / sqrt(1 - pow(eSpeed * sin(theta) / TMath::C(), 2)) - 1) *
              ME * TMath::C() * TMath::C()};

  const double z0_2{0.15};
  const double rho0_2{0.18};
  const double v0_2{-100};
  auto trap2 = new IdealPenningTrap(BFieldMag, v0_2, rho0_2, z0_2);
  const double C{2 * v0_2 / (z0_2 * z0_2 + 0.5 * rho0_2 * rho0_2)};
  const double z90{0.1};
  const double Ep90{-C * z90 * z90 * TMath::Qe() / 2};
  const double E90{Ep90 + ke90};
  const double Ep88{E90 - ke88};
  const double z88{sqrt(2 * Ep88 / (-TMath::Qe() * C))};
  std::cout << "90 degree potential energy = " << Ep90 / TMath::Qe() << " eV\n";
  std::cout << "88 degree potential energy = " << Ep88 / TMath::Qe() << " eV\n";
  std::cout << "z_88 = " << z88 * 100 << " cm\n";

  TVector3 startPos88(0, -rg, z88);
  TVector3 startPos90(0, -rg, z90);

  const double keDiff{ke88 - ke90};
  std::cout << "KE diff = " << keDiff / TMath::Qe() << " eV\n";

  TString trackFile88{outputDir + "/trackFile88.root"};
  ElectronTrajectoryGen traj88(trackFile88, trap2, startPos88, startVel88,
                               1e-12, 1e-6);
  // Open the input file and make various graphs
  auto fin88 = std::make_unique<TFile>(trackFile88, "read");
  auto tr88 = (TTree*)fin88->Get("tree");
  double time88{0};
  double xPos88{0}, yPos88{0}, zPos88{0};
  double xVel88{0}, yVel88{0}, zVel88{0};

  tr88->SetBranchAddress("time", &time88);
  tr88->SetBranchAddress("xPos", &xPos88);
  tr88->SetBranchAddress("yPos", &yPos88);
  tr88->SetBranchAddress("zPos", &zPos88);
  tr88->SetBranchAddress("xVel", &xVel88);
  tr88->SetBranchAddress("yVel", &yVel88);
  tr88->SetBranchAddress("zVel", &zVel88);
  auto grZ88 = new TGraph();
  setGraphAttr(grZ88);
  grZ88->SetLineColor(kBlue);
  grZ88->SetTitle("88 degrees; Time [s]; z [cm]");
  auto grX88 = new TGraph();
  setGraphAttr(grX88);
  grX88->SetLineColor(kBlue);
  grX88->SetTitle("88 degrees; Time [s]; x [cm]");

  auto grVr88 = new TGraph();
  setGraphAttr(grVr88);
  grVr88->SetLineColor(kBlue);
  grVr88->SetTitle("88 degrees; Time [s]; #beta_{r}");
  auto grVz88 = new TGraph();
  setGraphAttr(grVz88);
  grVz88->SetLineColor(kBlue);
  grVz88->SetTitle("88 degrees; Time [s]; #beta_{z}");

  auto grKE88 = new TGraph();
  setGraphAttr(grKE88);
  grKE88->SetLineColor(kBlue);
  grKE88->SetTitle("88 degrees; Time [s]; KE [eV]");

  for (int i{0}; i < tr88->GetEntries(); i++) {
    tr88->GetEntry(i);
    const double timeShift{27e-9 - 1.2e-11 + 17e-9};
    grX88->SetPoint(i, time88 + timeShift, xPos88 * 100);
    grZ88->SetPoint(i, time88 + timeShift, zPos88 * 100);
    double beta{sqrt(xVel88 * xVel88 + yVel88 * yVel88 + zVel88 * zVel88) /
                TMath::C()};
    double T{(1 / sqrt(1 - beta * beta) - 1) * ME * pow(TMath::C(), 2)};
    grKE88->SetPoint(i, time88 + timeShift, T / TMath::Qe());

    grVr88->SetPoint(i, time88 + timeShift,
                     sqrt(xVel88 * xVel88 + yVel88 * yVel88) / TMath::C());
    grVz88->SetPoint(i, time88 + timeShift, zVel88 / TMath::C());
  }
  delete tr88;
  fin88->Close();

  TString trackFile90{outputDir + "/trackFile90.root"};
  ElectronTrajectoryGen traj90(trackFile90, trap2, startPos90, startVel90,
                               1e-12, 1e-6);
  // Open the input file and make various graphs
  auto fin90 = std::make_unique<TFile>(trackFile90, "read");
  auto tr90 = (TTree*)fin90->Get("tree");
  double time90{0};
  double xPos90{0}, yPos90{0}, zPos90{0};
  double xVel90{0}, yVel90{0}, zVel90{0};
  tr90->SetBranchAddress("time", &time90);
  tr90->SetBranchAddress("xPos", &xPos90);
  tr90->SetBranchAddress("yPos", &yPos90);
  tr90->SetBranchAddress("zPos", &zPos90);
  tr90->SetBranchAddress("xVel", &xVel90);
  tr90->SetBranchAddress("yVel", &yVel90);
  tr90->SetBranchAddress("zVel", &zVel90);

  auto grZ90 = new TGraph();
  setGraphAttr(grZ90);
  grZ90->SetLineColor(kRed);
  grZ90->SetTitle("90 degrees; Time [s]; z [cm]");
  auto grX90 = new TGraph();
  setGraphAttr(grX90);
  grX90->SetLineColor(kRed);
  grX90->SetTitle("90 degrees; Time [s]; x [cm]");

  auto grVr90 = new TGraph();
  setGraphAttr(grVr90);
  grVr90->SetLineColor(kRed);
  grVr90->SetTitle("90 degrees; Time [s]; #beta_{r}");
  auto grVz90 = new TGraph();
  setGraphAttr(grVz90);
  grVz90->SetLineColor(kRed);
  grVz90->SetTitle("90 degrees; Time [s]; #beta_{z}");

  auto grKE90 = new TGraph();
  setGraphAttr(grKE90);
  grKE90->SetLineColor(kRed);
  grKE90->SetTitle("90 degrees; Time [s]; KE [eV]");

  for (int i{0}; i < tr90->GetEntries(); i++) {
    tr90->GetEntry(i);
    grX90->SetPoint(i, time90, xPos90 * 100);
    grZ90->SetPoint(i, time90, zPos90 * 100);
    double beta{sqrt(xVel90 * xVel90 + yVel90 * yVel90 + zVel90 * zVel90) /
                TMath::C()};
    double T{(1 / sqrt(1 - beta * beta) - 1) * ME * pow(TMath::C(), 2)};
    grKE90->SetPoint(i, time90, T / TMath::Qe());

    grVr90->SetPoint(i, time90,
                     sqrt(xVel90 * xVel90 + yVel90 * yVel90) / TMath::C());
    grVz90->SetPoint(i, time90, zVel90 / TMath::C());
  }
  delete tr90;
  fin90->Close();

  fout->cd();
  grX88->Write("grX88");
  grX90->Write("grX90");
  grZ88->Write("grZ88");
  grZ90->Write("grZ90");
  grKE88->Write("grKE88");
  grKE90->Write("grKE90");
  grVr88->Write("grVr88");
  grVr90->Write("grVr90");
  grVz88->Write("grVz88");
  grVz90->Write("grVz90");

  fout->Close();
  return 0;
}