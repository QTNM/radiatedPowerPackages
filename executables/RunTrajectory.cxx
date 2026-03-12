/*
  RunTrajectory.cxx

  Config-driven electron trajectory generation.
  Reads a YAML config file and generates a trajectory ROOT file.

  Usage: RunTrajectory <config.yaml>
*/

#include <cmath>
#include <ctime>
#include <iostream>
#include <string>
#include <tuple>

#include "TFile.h"
#include "TTree.h"
#include "TVector3.h"

#include "physics/ElectronDynamics/BorisSolver.h"
#include "utilities/BasicCore/Constants.h"
#include "utilities/ConfigParser/ConfigParser.h"

using namespace rad;

int main(int argc, char* argv[]) {
  if (argc != 2) {
    std::cerr << "Usage: " << argv[0] << " <config.yaml>" << std::endl;
    return 1;
  }

  try {
    auto config = config::LoadTrajectoryConfig(argv[1]);

    std::cout << "Output file: " << config.simulation.outputFile << std::endl;
    std::cout << "Simulation time: " << config.simulation.time << " s"
              << std::endl;
    std::cout << "Step size: " << config.simulation.stepSize << " s"
              << std::endl;
    std::cout << "Energy loss: "
              << (config.simulation.energyLoss ? "on" : "off") << std::endl;

    // Set up the Boris solver
    double tau = config.simulation.energyLoss ? 2 * R_E / (3 * C) : 0.0;
    BorisSolver solver(config.field.get(), -QE, ME, tau);

    // Open the output ROOT file
    TFile* fout = new TFile(config.simulation.outputFile.c_str(), "RECREATE");
    if (!fout || fout->IsZombie()) {
      throw std::runtime_error("Cannot create output file: " +
                               config.simulation.outputFile);
    }
    TTree* tree = new TTree("tree", "tree");

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
    TVector3 ePos = config.electron.position;
    TVector3 eVel = config.electron.velocity;
    TVector3 eAcc = solver.acc(ePos, eVel);

    time = config.simulation.initialTime;
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

    // Advance through the time steps
    const double stepSize = config.simulation.stepSize;
    const int nTimeSteps =
        static_cast<int>(std::round(config.simulation.time / stepSize));
    const auto& bounds = config.simulation.bounds;

    const clock_t beginTime = clock();
    const double printoutInterval = 1e-6;
    double printoutTime = printoutInterval;

    for (int i = 1; i < nTimeSteps; i++) {
      time = config.simulation.initialTime + double(i) * stepSize;
      std::tuple<TVector3, TVector3> outputStep =
          solver.advance_step(stepSize, ePos, eVel);

      if (time >= printoutTime) {
        std::cout << printoutTime << " seconds of trajectory simulated..."
                  << std::endl;
        printoutTime += printoutInterval;
      }

      ePos = std::get<0>(outputStep);
      eVel = std::get<1>(outputStep);
      eAcc = solver.acc(ePos, eVel);

      // Check geometric bounds
      double r = std::sqrt(ePos.X() * ePos.X() + ePos.Y() * ePos.Y());
      if (ePos.Z() < bounds.zMin || ePos.Z() > bounds.zMax ||
          r > bounds.rMax) {
        std::cout << "Electron exited bounds at t=" << time << " s"
                  << std::endl;
        break;
      }

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
    tree->Write("", TObject::kOverwrite);
    fout->Close();
    delete fout;

    const clock_t endTime = clock();
    std::cout << "Execution time: "
              << float(endTime - beginTime) / CLOCKS_PER_SEC << " seconds"
              << std::endl;

  } catch (const std::exception& e) {
    std::cerr << "Error: " << e.what() << std::endl;
    return 1;
  }

  return 0;
}
