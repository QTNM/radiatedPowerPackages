/*
  RunTrajectory.cxx

  Config-driven electron trajectory generation.
  Reads a YAML config file and generates trajectory ROOT file(s).
  Supports single or multiple electrons per config.

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

/// Replace all occurrences of "{index}" in str with the given index value
static std::string SubstituteIndex(const std::string& str, size_t index) {
  std::string result = str;
  const std::string placeholder = "{index}";
  size_t pos = 0;
  while ((pos = result.find(placeholder, pos)) != std::string::npos) {
    result.replace(pos, placeholder.length(), std::to_string(index));
    pos += std::to_string(index).length();
  }
  return result;
}

/// Run a single trajectory simulation and write to a ROOT file
static void RunSingleTrajectory(const std::string& outputFile,
                                BorisSolver& solver,
                                const config::ElectronConfig& electron,
                                const config::SimulationConfig& simulation) {
  TFile* fout = new TFile(outputFile.c_str(), "RECREATE");
  if (!fout || fout->IsZombie()) {
    throw std::runtime_error("Cannot create output file: " + outputFile);
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
  TVector3 ePos = electron.position;
  TVector3 eVel = electron.velocity;
  TVector3 eAcc = solver.acc(ePos, eVel);

  time = simulation.initialTime;
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
  const double stepSize = simulation.stepSize;
  const int nTimeSteps =
      static_cast<int>(std::round(simulation.time / stepSize));
  const auto& bounds = simulation.bounds;

  const double printoutInterval = 1e-6;
  double printoutTime = printoutInterval;

  for (int i = 1; i < nTimeSteps; i++) {
    time = simulation.initialTime + double(i) * stepSize;
    std::tuple<TVector3, TVector3> outputStep =
        solver.advance_step(stepSize, ePos, eVel);

    if (time >= printoutTime) {
      std::cout << "  " << printoutTime
                << " seconds of trajectory simulated..." << std::endl;
      printoutTime += printoutInterval;
    }

    ePos = std::get<0>(outputStep);
    eVel = std::get<1>(outputStep);
    eAcc = solver.acc(ePos, eVel);

    // Check geometric bounds
    double r = std::sqrt(ePos.X() * ePos.X() + ePos.Y() * ePos.Y());
    if (ePos.Z() < bounds.zMin || ePos.Z() > bounds.zMax || r > bounds.rMax) {
      std::cout << "  Electron exited bounds at t=" << time << " s"
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
}

int main(int argc, char* argv[]) {
  if (argc != 2) {
    std::cerr << "Usage: " << argv[0] << " <config.yaml>" << std::endl;
    return 1;
  }

  try {
    auto config = config::LoadTrajectoryConfig(argv[1]);

    std::cout << "Simulation time: " << config.simulation.time << " s"
              << std::endl;
    std::cout << "Step size: " << config.simulation.stepSize << " s"
              << std::endl;
    std::cout << "Energy loss: "
              << (config.simulation.energyLoss ? "on" : "off") << std::endl;
    std::cout << "Electrons: " << config.electrons.size() << std::endl;

    // Validate output filename for multi-electron runs
    if (config.electrons.size() > 1 &&
        config.simulation.outputFile.find("{index}") == std::string::npos) {
      throw std::runtime_error(
          "Multi-electron config requires '{index}' placeholder in output_file "
          "(e.g. \"traj_{index}.root\")");
    }

    // Set up the Boris solver (shared across all electrons — same field)
    double tau = config.simulation.energyLoss ? 2 * R_E / (3 * C) : 0.0;
    BorisSolver solver(config.field.get(), -QE, ME, tau);

    const clock_t beginTime = clock();

    for (size_t i = 0; i < config.electrons.size(); i++) {
      std::string filename =
          SubstituteIndex(config.simulation.outputFile, i);

      if (config.electrons.size() > 1) {
        std::cout << "Electron " << i + 1 << "/" << config.electrons.size()
                  << ": " << filename << std::endl;
      } else {
        std::cout << "Output file: " << filename << std::endl;
      }

      RunSingleTrajectory(filename, solver, config.electrons[i],
                           config.simulation);
    }

    const clock_t endTime = clock();
    std::cout << "Total execution time: "
              << float(endTime - beginTime) / CLOCKS_PER_SEC << " seconds"
              << std::endl;

  } catch (const std::exception& e) {
    std::cerr << "Error: " << e.what() << std::endl;
    return 1;
  }

  return 0;
}
