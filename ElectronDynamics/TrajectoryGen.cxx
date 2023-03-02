/// TrajectoryGen.cxx

#include "ElectronDynamics/TrajectoryGen.h"

#include <cmath>
#include <memory>
#include <tuple>

#include "ElectronDynamics/BaseField.h"
#include "ElectronDynamics/BorisSolver.h"
#include "TFile.h"
#include "TMath.h"
#include "TString.h"
#include "TTree.h"
#include "TVector3.h"

rad::ElectronTrajectoryGen::ElectronTrajectoryGen(
    TString outputFile, BaseField *field, TVector3 initPos, TVector3 initVel,
    double simStepSize, double simTime, double initialSimTime, double tau) {
  // Check the file path can be opened in
  auto foutTest = std::make_unique<TFile>(outputFile, "RECREATE");
  if (!foutTest) {
    // File path not opened correctly
    std::cout << "File cannot be created. Exiting..." << std::endl;
    exit(1);
  } else {
    // File path is all good
    foutTest->Close();
  }

  solver = BorisSolver(field, -TMath::Qe(), ME, tau);

  // Check that various input values make sense
  if (simStepSize <= 0) {
    std::cout << "Invalid simulation step size (" << simStepSize
              << "). Exiting..." << std::endl;
    exit(1);
  }
  if (simTime <= 0) {
    std::cout << "Invalid simulation time (" << simTime << "). Exiting..."
              << std::endl;
    exit(1);
  }

  // Open the output file
  auto fout = std::make_unique<TFile>(outputFile, "RECREATE");
  auto tree = std::make_unique<TTree>("tree", "tree");
  double time{};
  double xPos{}, yPos{}, zPos{};
  double xVel{}, yVel{}, zVel{};
  double xAcc{}, yAcc{}, zAcc{};
  tree->Branch("time", &time);
  tree->Branch("xPos", &xPos);
  tree->Branch("yPos", &yPos);
  tree->Branch("zPos", &zPos);
  tree->Branch("xVel", &xVel);
  tree->Branch("yVel", &yVel);
  tree->Branch("zVel", &zVel);
  tree->Branch("xAcc", &xAcc);
  tree->Branch("yAcc", &yAcc);
  tree->Branch("zAcc", &zAcc);

  // Set the initial state
  TVector3 eAcc = solver.acc(initPos, initVel);
  time = initialSimTime;
  xPos = initPos.X();
  yPos = initPos.Y();
  zPos = initPos.Z();
  xVel = initPos.X();
  yVel = initPos.Y();
  zVel = initPos.Z();
  xAcc = eAcc.X();
  yAcc = eAcc.Y();
  zAcc = eAcc.Z();
  tree->Fill();

  TVector3 ePos = initPos;
  TVector3 eVel = initPos;

  double nTimeSteps{round(simTime / simStepSize)};
  // Advance through the time steps
  const double printoutTime{1e-6};  // seconds
  for (int i = 1; i < nTimeSteps; i++) {
    time = initialSimTime + double(i) * simStepSize;
    std::tuple<TVector3, TVector3> outputStep =
        solver.advance_step(simStepSize, ePos, eVel);

    if (std::fmod(time, printoutTime) < simStepSize) {
      std::cout << time << " seconds of trajectory simulated..." << std::endl;
    }

    ePos = std::get<0>(outputStep);
    eVel = std::get<1>(outputStep);
    eAcc = solver.acc(ePos, eVel);

    xPos = ePos.X();
    yPos = ePos.Y();
    zPos = ePos.Z();
    xVel = eVel.X();
    yVel = eVel.Y();
    zVel = eVel.Z();
    xAcc = eAcc.X();
    yAcc = eAcc.Y();
    zAcc = eAcc.Z();

    tree->Fill();
  }
  fout->cd();
  tree->Write();
  fout->Close();
}