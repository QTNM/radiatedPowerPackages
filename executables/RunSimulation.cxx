/*
  RunSimulation.cxx

  Config-driven end-to-end signal generation pipeline.
  Reads a YAML pipeline config, generates electron trajectories, computes
  downmixed signals via the configured detector, and writes outputs to HDF5.

  Usage: RunSimulation <config.yaml>
*/

#include <boost/lexical_cast.hpp>
#include <boost/uuid/uuid.hpp>
#include <boost/uuid/uuid_generators.hpp>
#include <boost/uuid/uuid_io.hpp>
#include <cmath>
#include <ctime>
#include <filesystem>
#include <iostream>
#include <string>
#include <tuple>
#include <variant>
#include <vector>

#include "H5Cpp.h"
#include "TFile.h"
#include "TGraph.h"
#include "TTree.h"
#include "TVector3.h"

#include "physics/ElectronDynamics/BorisSolver.h"
#include "physics/SignalProcessing/LocalOscillator.h"
#include "physics/SignalProcessing/NoiseFunc.h"
#include "physics/SignalProcessing/Signal.h"
#include "utilities/BasicCore/Constants.h"
#include "utilities/ConfigParser/PipelineConfig.h"
#include "utilities/HDF5Writer/HDF5Writer.h"

using namespace rad;

static std::string MakeUUID() {
  return boost::lexical_cast<std::string>(
      (boost::uuids::random_generator())());
}

/// Run a single trajectory and write result to a ROOT file
static void RunSingleTrajectory(const std::string& outputFile,
                                BorisSolver& solver,
                                const config::ElectronConfig& electron,
                                const config::SimulationConfig& simulation) {
  TFile* fout = new TFile(outputFile.c_str(), "RECREATE");
  if (!fout || fout->IsZombie()) {
    throw std::runtime_error("Cannot create trajectory file: " + outputFile);
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

  TVector3 ePos = electron.position;
  TVector3 eVel = electron.velocity;
  TVector3 eAcc = solver.acc(ePos, eVel);

  time = simulation.initialTime;
  xPos = ePos.X(); yPos = ePos.Y(); zPos = ePos.Z();
  xVel = eVel.X(); yVel = eVel.Y(); zVel = eVel.Z();
  xAcc = eAcc.X(); yAcc = eAcc.Y(); zAcc = eAcc.Z();
  tree->Fill();

  const double stepSize = simulation.stepSize;
  const int nSteps = static_cast<int>(std::round(simulation.time / stepSize));
  const auto& bounds = simulation.bounds;

  const double printInterval = 1e-6;
  double printTime = printInterval;

  for (int i = 1; i < nSteps; i++) {
    time = simulation.initialTime + double(i) * stepSize;
    std::tuple<TVector3, TVector3> step =
        solver.advance_step(stepSize, ePos, eVel);

    if (time >= printTime) {
      std::cout << "  " << printTime << " s simulated...\n";
      printTime += printInterval;
    }

    ePos = std::get<0>(step);
    eVel = std::get<1>(step);
    eAcc = solver.acc(ePos, eVel);

    double r = std::sqrt(ePos.X() * ePos.X() + ePos.Y() * ePos.Y());
    if (ePos.Z() < bounds.zMin || ePos.Z() > bounds.zMax || r > bounds.rMax) {
      std::cout << "  Electron exited bounds at t=" << time << " s\n";
      break;
    }

    xPos = ePos.X(); yPos = ePos.Y(); zPos = ePos.Z();
    xVel = eVel.X(); yVel = eVel.Y(); zVel = eVel.Z();
    xAcc = eAcc.X(); yAcc = eAcc.Y(); zAcc = eAcc.Z();
    tree->Fill();
  }

  fout->cd();
  tree->Write("", TObject::kOverwrite);
  fout->Close();
  delete fout;
}

/// Build GaussianNoise terms from config
static std::vector<GaussianNoise> BuildNoiseTerms(
    const std::vector<config::NoiseConfig>& noiseCfg) {
  std::vector<GaussianNoise> noise;
  for (const auto& nc : noiseCfg) {
    noise.emplace_back(nc.temperature, nc.resistance, nc.seed);
  }
  return noise;
}

/// Write requested signal datasets into an HDF5 group
static void WriteSignalOutputs(H5::Group& group, Signal& sig,
                                const std::vector<std::string>& datasets) {
  for (const auto& dsName : datasets) {
    if (dsName == "voltage_i") {
      hdf5::WriteGraphDataset(group, "voltage_i", sig.GetVITimeDomain());
    } else if (dsName == "voltage_q") {
      hdf5::WriteGraphDataset(group, "voltage_q", sig.GetVQTimeDomain());
    } else if (dsName == "power_spectrum") {
      constexpr double kLoadOhm = 50.0;
      TGraph* ps = sig.GetVIPowerPeriodogram(kLoadOhm);
      hdf5::WriteGraphDataset(group, "power_spectrum", ps,
                              /*writeTimeStep=*/false);
      delete ps;
    } else {
      std::cerr << "  Warning: unknown dataset '" << dsName
                << "' — skipping\n";
    }
  }
}

/// Write simulation metadata as HDF5 attributes on a group
static void WriteMetadata(H5::Group& group,
                          const config::SimulationConfig& sim,
                          const config::SignalProcessingConfig& sp) {
  hdf5::WriteScalarAttribute(group, "simulation_time_s", sim.time);
  hdf5::WriteScalarAttribute(group, "step_size_s", sim.stepSize);
  hdf5::WriteScalarAttribute(group, "sample_rate_hz", sp.sampleRate);
  hdf5::WriteScalarAttribute(group, "lo_frequency_hz", sp.loFrequency);
  if (sp.acquisitionTime > 0) {
    hdf5::WriteScalarAttribute(group, "acquisition_time_s",
                               sp.acquisitionTime);
  }
}

/// Process one electron: generate trajectory, compute signal(s), write HDF5
static void ProcessElectron(H5::H5File& h5file, const std::string& groupName,
                             BorisSolver& solver,
                             const config::ElectronConfig& electron,
                             const config::PipelineConfig& config,
                             const LocalOscillator& lo,
                             const std::vector<GaussianNoise>& noiseTerms,
                             const std::string& trajFile) {
  RunSingleTrajectory(trajFile, solver, electron, config.simulation);

  TString trajPath(trajFile.c_str());
  const auto& sp = config.signalProcessing;
  const auto& output = config.output;

  H5::Group group(h5file.createGroup(groupName));

  std::visit(
      [&](const auto& det) {
        using T = std::decay_t<decltype(det)>;

        if constexpr (std::is_same_v<T, config::AntennaDetector>) {
          // Build raw pointer vector for Signal
          std::vector<IAntenna*> antPtrs;
          for (const auto& a : det.antennas) {
            antPtrs.push_back(a.get());
          }
          Signal sig(trajPath, antPtrs, lo, sp.sampleRate, noiseTerms,
                     sp.acquisitionTime);
          WriteSignalOutputs(group, sig, output.datasets);

        } else if constexpr (std::is_same_v<T, config::WaveguideDetector>) {
          for (size_t j = 0; j < det.probes.size(); j++) {
            Signal sig(trajPath, det.waveguide.get(), lo, sp.sampleRate,
                       det.probes[j], noiseTerms, sp.acquisitionTime);

            if (det.probes.size() == 1) {
              WriteSignalOutputs(group, sig, output.datasets);
            } else {
              std::string probeName = "probe_" + std::to_string(j);
              H5::Group probeGroup(group.createGroup(probeName));
              WriteSignalOutputs(probeGroup, sig, output.datasets);
            }
          }

        } else if constexpr (std::is_same_v<T, config::CavityDetector>) {
          Signal sig(trajPath, det.cavity.get(), lo, sp.sampleRate, noiseTerms,
                     sp.acquisitionTime);
          WriteSignalOutputs(group, sig, output.datasets);
        }
      },
      config.detector);

  if (output.metadata) {
    WriteMetadata(group, config.simulation, sp);
  }
}

int main(int argc, char* argv[]) {
  if (argc != 2) {
    std::cerr << "Usage: " << argv[0] << " <config.yaml>\n";
    return 1;
  }

  try {
    auto config = config::LoadPipelineConfig(argv[1]);

    const auto& sim = config.simulation;
    const auto& sp = config.signalProcessing;
    const auto& output = config.output;

    std::cout << "Simulation time:  " << sim.time << " s\n";
    std::cout << "Step size:        " << sim.stepSize << " s\n";
    std::cout << "Energy loss:      " << (sim.energyLoss ? "on" : "off")
              << "\n";
    std::cout << "Electrons:        " << config.electrons.size() << "\n";
    std::cout << "Sample rate:      " << sp.sampleRate / 1e9 << " GHz\n";
    std::cout << "LO frequency:     " << sp.loFrequency / 1e9 << " GHz\n";
    std::cout << "Output file:      " << output.file << "\n";

    // Set up shared physics objects
    double tau = sim.energyLoss ? 2.0 * R_E / (3.0 * C) : 0.0;
    BorisSolver solver(config.field.get(), -QE, ME, tau);

    LocalOscillator lo(2.0 * PI * sp.loFrequency);
    std::vector<GaussianNoise> noiseTerms = BuildNoiseTerms(sp.noise);

    // Create temp trajectory directory
    std::string trajDir = output.trajectoryDir + "/RunSimulation_" + MakeUUID();
    std::filesystem::create_directories(trajDir);

    // Open HDF5 output file
    H5::H5File h5file(output.file, H5F_ACC_TRUNC);

    const clock_t startTime = clock();
    const size_t nElectrons = config.electrons.size();
    const bool multiElectron = nElectrons > 1;

    for (size_t i = 0; i < nElectrons; i++) {
      std::string trajFile =
          trajDir + "/traj_" + std::to_string(i) + ".root";
      std::string groupName =
          multiElectron ? output.group + "/electron_" + std::to_string(i)
                        : output.group;

      if (multiElectron) {
        std::cout << "\nElectron " << i + 1 << "/" << nElectrons << "\n";
      }

      ProcessElectron(h5file, groupName, solver, config.electrons[i], config,
                      lo, noiseTerms, trajFile);

      if (output.cleanupTrajectories) {
        std::filesystem::remove(trajFile);
      }
    }

    // Clean up temp directory (empty if trajectories were removed)
    if (output.cleanupTrajectories) {
      std::filesystem::remove(trajDir);
    }

    const clock_t endTime = clock();
    std::cout << "\nTotal time: "
              << float(endTime - startTime) / CLOCKS_PER_SEC << " s\n";

  } catch (const std::exception& e) {
    std::cerr << "Error: " << e.what() << "\n";
    return 1;
  }

  return 0;
}
