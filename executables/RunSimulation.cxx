/*
  RunSimulation.cxx

  Config-driven end-to-end signal generation pipeline.
  Reads a YAML pipeline config, generates electron trajectories (optionally
  with gas scattering), computes downmixed signals via the configured detector,
  and writes outputs to HDF5.

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
#include <memory>
#include <random>
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
#include "physics/Scattering/ElasticScatter.h"
#include "physics/Scattering/InelasticScatter.h"
#include "physics/SignalProcessing/LocalOscillator.h"
#include "physics/SignalProcessing/NoiseFunc.h"
#include "physics/SignalProcessing/Signal.h"
#include "utilities/BasicCore/Constants.h"
#include "utilities/BasicCore/Physics.h"
#include "utilities/ConfigParser/PipelineConfig.h"
#include "utilities/HDF5Writer/HDF5Writer.h"

using namespace rad;

static std::string MakeUUID() {
  return boost::lexical_cast<std::string>(
      (boost::uuids::random_generator())());
}

// ─── Species mapping ─────────────────────────────────────────────────────────

struct SpeciesPhysics {
  Species inelasticSpecies;
  unsigned int atomicNumber;
  unsigned int atomicMass;
};

/// Map a species name to its elastic/inelastic physics parameters.
/// T/T2 use the same inelastic cross-sections as H/H2 respectively.
static SpeciesPhysics MapSpecies(const std::string& name) {
  if (name == "H")  return {H,  1, 1};
  if (name == "H2") return {H2, 2, 2};
  if (name == "He") return {He, 2, 4};
  if (name == "T")  return {H,  1, 3};
  if (name == "T2") return {H2, 2, 6};
  throw std::runtime_error("Unknown gas species: " + name);
}

// ─── Trajectory generation ───────────────────────────────────────────────────

/// Run a single trajectory and write to a ROOT file.
/// When scattering is present the loop is event-driven; otherwise it
/// runs fixed Boris steps — both cases produce the same ROOT output format.
static void RunSingleTrajectory(
    const std::string& outputFile, BorisSolver& solver,
    const config::ElectronConfig& electron,
    const config::SimulationConfig& simulation,
    const std::optional<config::ScatteringConfig>& scattering) {

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

  auto Record = [&](double t, const TVector3& p, const TVector3& v,
                    const TVector3& a) {
    time = t;
    xPos = p.X(); yPos = p.Y(); zPos = p.Z();
    xVel = v.X(); yVel = v.Y(); zVel = v.Z();
    xAcc = a.X(); yAcc = a.Y(); zAcc = a.Z();
    tree->Fill();
  };

  TVector3 ePos = electron.position;
  TVector3 eVel = electron.velocity;
  TVector3 eAcc = solver.acc(ePos, eVel);
  Record(simulation.initialTime, ePos, eVel, eAcc);

  const double stepSize = simulation.stepSize;
  const double maxTime  = simulation.initialTime + simulation.time;
  const auto&  bounds   = simulation.bounds;

  const double printInterval = 1e-6;
  double printTime = simulation.initialTime + printInterval;

  std::mt19937 gen(std::random_device{}());

  double t = simulation.initialTime;

  while (t < maxTime) {
    double dt       = stepSize;
    bool scattered  = false;

    // ── Per-species scattering objects (only allocated when needed) ─────────
    std::vector<std::unique_ptr<ElasticScatter>>   elastics;
    std::vector<std::unique_ptr<InelasticScatter>> inelastics;
    std::vector<double> elasticXSecs, inelasticXSecs, rates;
    double totalRate = 0.0;

    if (scattering) {
      double gamma = 1.0 / std::sqrt(1.0 - eVel.Mag2() / (C * C));
      double ke    = (gamma - 1.0) * ME_EV;
      double vMag  = eVel.Mag();

      for (const auto& gas : scattering->gases) {
        auto sp = MapSpecies(gas.species);
        auto el = std::make_unique<ElasticScatter>(ke, sp.atomicNumber,
                                                   sp.atomicMass);
        auto in = std::make_unique<InelasticScatter>(ke, sp.inelasticSpecies);

        double eSec = el->GetTotalXSec();
        double iSec = in->GetTotalXSec();
        double rate = gas.density * (eSec + iSec) * vMag;

        elasticXSecs.push_back(eSec);
        inelasticXSecs.push_back(iSec);
        rates.push_back(rate);
        totalRate += rate;

        elastics.push_back(std::move(el));
        inelastics.push_back(std::move(in));
      }

      std::exponential_distribution<double> scatterDist(totalRate);
      double scatterTime = scatterDist(gen);

      if (scatterTime < stepSize) {
        dt        = scatterTime;
        scattered = true;
      }
    }

    // ── Advance particle ────────────────────────────────────────────────────
    auto step = solver.advance_step(dt, ePos, eVel);
    ePos = std::get<0>(step);
    eVel = std::get<1>(step);
    t += dt;

    double r = std::sqrt(ePos.X() * ePos.X() + ePos.Y() * ePos.Y());
    bool outOfBounds = (ePos.Z() < bounds.zMin || ePos.Z() > bounds.zMax ||
                        r > bounds.rMax);

    eAcc = solver.acc(ePos, eVel);
    Record(t, ePos, eVel, eAcc);

    if (outOfBounds) {
      std::cout << "  Electron exited bounds at t=" << t << " s\n";
      break;
    }

    // ── Apply scatter kinematics ─────────────────────────────────────────────
    if (scattered) {
      std::uniform_real_distribution<double> uni(0.0, 1.0);

      // Pick which species
      double rSpecies = uni(gen) * totalRate;
      size_t j = rates.size() - 1;
      double cumRate = 0.0;
      for (size_t k = 0; k < rates.size(); k++) {
        cumRate += rates[k];
        if (rSpecies <= cumRate) { j = k; break; }
      }

      double scatterAngle = 0.0;
      if (uni(gen) < elasticXSecs[j] / (elasticXSecs[j] + inelasticXSecs[j])) {
        // Elastic
        scatterAngle     = elastics[j]->GetRandomScatteringAngle();
        double newKE     = elastics[j]->GetEnergyAfterScatter(scatterAngle);
        eVel = elastics[j]->GetScatteredVector(eVel, newKE, scatterAngle);
      } else {
        // Inelastic
        double W     = inelastics[j]->GetRandomW();
        scatterAngle = inelastics[j]->GetRandomTheta(W);
        eVel = inelastics[j]->GetScatteredVector(eVel, W, scatterAngle);
      }

      double pitchDeg =
          std::abs(std::atan(eVel.Perp() / eVel.Z())) * 180.0 / M_PI;
      std::cout << "  Scatter at t=" << t << " s"
                << "  angle=" << scatterAngle * 180.0 / M_PI << " deg"
                << "  pitch=" << pitchDeg << " deg\n";

      eAcc = solver.acc(ePos, eVel);
      Record(t, ePos, eVel, eAcc);
    }

    if (t >= printTime) {
      std::cout << "  " << printTime << " s simulated...\n";
      printTime += printInterval;
    }
  }

  fout->cd();
  tree->Write("", TObject::kOverwrite);
  fout->Close();
  delete fout;
}

// ─── Signal construction helpers ─────────────────────────────────────────────

static std::vector<GaussianNoise> BuildNoiseTerms(
    const std::vector<config::NoiseConfig>& noiseCfg) {
  std::vector<GaussianNoise> noise;
  for (const auto& nc : noiseCfg) {
    noise.emplace_back(nc.temperature, nc.resistance, nc.seed);
  }
  return noise;
}

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
      std::cerr << "  Warning: unknown dataset '" << dsName << "' — skipping\n";
    }
  }
}

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

// ─── Per-electron pipeline ────────────────────────────────────────────────────

static void ProcessElectron(H5::H5File& h5file, const std::string& groupName,
                             BorisSolver& solver,
                             const config::ElectronConfig& electron,
                             const config::PipelineConfig& config,
                             const LocalOscillator& lo,
                             const std::vector<GaussianNoise>& noiseTerms,
                             const std::string& trajFile) {
  RunSingleTrajectory(trajFile, solver, electron, config.simulation,
                      config.scattering);

  TString trajPath(trajFile.c_str());
  const auto& sp     = config.signalProcessing;
  const auto& output = config.output;

  H5::Group group(h5file.createGroup(groupName));

  std::visit(
      [&](const auto& det) {
        using T = std::decay_t<decltype(det)>;

        if constexpr (std::is_same_v<T, config::AntennaDetector>) {
          std::vector<IAntenna*> antPtrs;
          for (const auto& a : det.antennas) antPtrs.push_back(a.get());
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
              H5::Group probeGroup(
                  group.createGroup("probe_" + std::to_string(j)));
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

// ─── Main ─────────────────────────────────────────────────────────────────────

int main(int argc, char* argv[]) {
  if (argc != 2) {
    std::cerr << "Usage: " << argv[0] << " <config.yaml>\n";
    return 1;
  }

  try {
    auto config = config::LoadPipelineConfig(argv[1]);

    const auto& sim    = config.simulation;
    const auto& sp     = config.signalProcessing;
    const auto& output = config.output;

    std::cout << "Simulation time:  " << sim.time << " s\n";
    std::cout << "Step size:        " << sim.stepSize << " s\n";
    std::cout << "Energy loss:      " << (sim.energyLoss ? "on" : "off") << "\n";
    std::cout << "Electrons:        " << config.electrons.size() << "\n";
    std::cout << "Sample rate:      " << sp.sampleRate / 1e9 << " GHz\n";
    std::cout << "LO frequency:     " << sp.loFrequency / 1e9 << " GHz\n";
    std::cout << "Output file:      " << output.file << "\n";
    if (config.scattering) {
      std::cout << "Scattering:       on (" << config.scattering->gases.size()
                << " gas species)\n";
    } else {
      std::cout << "Scattering:       off\n";
    }

    double tau = sim.energyLoss ? 2.0 * R_E / (3.0 * C) : 0.0;
    BorisSolver solver(config.field.get(), -QE, ME, tau);

    LocalOscillator lo(2.0 * PI * sp.loFrequency);
    std::vector<GaussianNoise> noiseTerms = BuildNoiseTerms(sp.noise);

    std::string trajDir =
        output.trajectoryDir + "/RunSimulation_" + MakeUUID();
    std::filesystem::create_directories(trajDir);

    H5::H5File h5file(output.file, H5F_ACC_TRUNC);

    const clock_t startTime  = clock();
    const size_t  nElectrons = config.electrons.size();
    const bool    multi      = nElectrons > 1;

    for (size_t i = 0; i < nElectrons; i++) {
      std::string trajFile =
          trajDir + "/traj_" + std::to_string(i) + ".root";
      std::string groupName =
          multi ? output.group + "/electron_" + std::to_string(i)
                : output.group;

      if (multi) {
        std::cout << "\nElectron " << i + 1 << "/" << nElectrons << "\n";
      }

      ProcessElectron(h5file, groupName, solver, config.electrons[i], config,
                      lo, noiseTerms, trajFile);

      if (output.cleanupTrajectories) {
        std::filesystem::remove(trajFile);
      }
    }

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
