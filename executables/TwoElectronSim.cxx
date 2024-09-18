/*
  TwoElectronSim.cxx

  Test program to look at the motion of two electrons in a magnetic trap.
*/

#include <getopt.h>
#include <unistd.h>

#include <boost/lexical_cast.hpp>
#include <boost/uuid/uuid.hpp>
#include <boost/uuid/uuid_generators.hpp>
#include <boost/uuid/uuid_io.hpp>
#include <ctime>
#include <iostream>
#include <memory>
#include <random>
#include <vector>

#include "BasicFunctions/BasicFunctions.h"
#include "BasicFunctions/Constants.h"
#include "BasicFunctions/TritiumSpectrum.h"
#include "ElectronDynamics/BorisSolver.h"
#include "ElectronDynamics/ManyParticleSolver.h"
#include "ElectronDynamics/Particle.h"
#include "ElectronDynamics/QTNMFields.h"
#include "ElectronDynamics/TrajectoryGen.h"
#include "TFile.h"
#include "TGraph.h"
#include "TString.h"
#include "TTree.h"
#include "TVector3.h"

using namespace rad;

typedef std::piecewise_linear_distribution<double> pld;

/// @brief Function for generating a random file string
/// @return A unique string which can be used for file names
std::string make_uuid() {
  return boost::lexical_cast<std::string>((boost::uuids::random_generator())());
}

pld GenerateSpectrum(int NPoints) {
  std::vector<double> energyVec{};
  std::vector<double> decayRateVec{};
  const double endpointEnergy{18.6e3};

  for (int i{0}; i < NPoints; i++) {
    double E{0.1 + endpointEnergy * double(i) / double(NPoints - 1)};
    energyVec.push_back(E);
    decayRateVec.push_back(TritiumDecayRate(E, 0, 0, 0, endpointEnergy));
  }

  return pld(energyVec.begin(), energyVec.end(), decayRateVec.begin());
}

/// @brief Generate a random direction vector for
/// @param mt Random number generator
/// @return A TVector3 with a random direction
TVector3 GenerateRandomDirection(std::mt19937& mt) {
  std::uniform_real_distribution<double> uniDist(0, 1);
  double u1{uniDist(mt)};
  double u2{uniDist(mt)};
  double phiGen{u1 * 2 * M_PI};
  double thetaGen{acos(2 * u2 - 1)};
  return TVector3(sin(thetaGen) * cos(phiGen), sin(thetaGen) * sin(phiGen),
                  cos(thetaGen));
}

/// @brief Generates a random velocity from a inputted ke distribution
/// @param p PDF of kinetic energies
/// @param mt Random number generator
/// @return A velocity vector in units of m s^-1
TVector3 GenerateDecayVelocity(pld& p, std::mt19937& mt) {
  double ke{p(mt)};
  // Calculate electron speed
  double v{GetSpeedFromKE(ke, ME)};

  return v * GenerateRandomDirection(mt);
}

int main(int argc, char* argv[]) {
  // Parse the command line arguments
  std::string outputDir;
  bool drawFromSpectrum{false};
  int opt{};
  int nSims{1};
  while ((opt = getopt(argc, argv, ":o:n:r")) != -1) {
    switch (opt) {
      case 'o':
        std::cout << "Output directory: " << optarg << std::endl;
        outputDir = optarg;
        break;
      case 'r':
        std::cout << "Drawing from spectrum" << std::endl;
        drawFromSpectrum = true;
        break;
      case 'n':
        std::cout << "Number of simulations: " << optarg << std::endl;
        nSims = std::stoi(optarg);
        break;
      case ':':
        std::cerr << "Option -" << optopt << " requires an argument.\n";
        return 1;
      case '?':
        std::cerr << "Unknown option -" << optopt << std::endl;
        return 1;
    }
  }

  std::string outFileName{"/out_" + make_uuid() + ".root"};
  TFile outFile((outputDir + outFileName).c_str(), "RECREATE");

  // We want to look at the effect on an endpoint electron
  const double endpointEnergy{18.575e3};                           // eV
  const double endpointSpeed{GetSpeedFromKE(endpointEnergy, ME)};  // m/s

  // Magnetic field setup
  const double B0{1.0};                                                  // T
  const double trapDepth{4e-3};                                          // T
  const double trapVolume{1e-3};                                         // m^-3
  const double trapLength{0.15};                                         // m
  const double trapA{trapVolume / trapLength};                           // m^2
  const double trapCoilRadius{sqrt(trapVolume / M_PI)};                  // m
  const double trapCoilCurrent{2.0 * trapDepth * trapCoilRadius / MU0};  // A
  auto bathtubField{new BathtubField(trapCoilRadius, trapCoilCurrent,
                                     -trapLength / 2, trapLength / 2,
                                     TVector3(0, 0, B0))};

  // Simulation timings setup
  double cycFreq0eV{CalcCyclotronFreq(0, B0)};
  const double dt{1 / (12 * cycFreq0eV)};  // s
  const double tMax{10e-6};                // s

  // Calculate minimum trapping angle
  const double dThetaMax{acos(sqrt(B0 / (B0 + trapDepth)))};
  const double thetaMin{M_PI / 2 - dThetaMax};
  std::cout << "Minimum trapping angle: " << thetaMin * 180 / M_PI << std::endl;

  // Make a plot of the field
  const uint nPoints{300};
  const double zMin{-trapLength / 2 * 1.1};
  const double zMax{trapLength / 2 * 1.1};
  TGraph* fieldPlot = new TGraph();
  setGraphAttr(fieldPlot);
  fieldPlot->GetXaxis()->SetTitle("z [m]");
  fieldPlot->GetYaxis()->SetTitle("|B| [T]");
  for (uint i{0}; i < nPoints; i++) {
    double z{zMin + double(i) * (zMax - zMin) / double(nPoints - 1)};
    TVector3 pos(0, 0, z);
    fieldPlot->SetPoint(i, z, bathtubField->evaluate_field_magnitude(pos));
  }
  outFile.cd();
  fieldPlot->Write("fieldPlot");

  pld p{GenerateSpectrum(1000)};

  int simCounter{0};

  // Create a TTree for writing the results of the sims
  TTree simTree("simTree", "Simulation results");
  double endpointStartPos[3];
  double endpointStartVel[3];
  double endpointKEInit;
  double endpointKE_u;
  double endpointKE_p;
  double endpointPitchAngle;
  double perturberStartPos[3];
  double perturberStartVel[3];
  double perturberKEInit;
  double perturberPitchAngle;
  simTree.Branch("endpointStartPos", endpointStartPos, "endpointStartPos[3]/D");
  simTree.Branch("endpointStartVel", endpointStartVel, "endpointStartVel[3]/D");
  simTree.Branch("endpointKEInit", &endpointKEInit, "endpointKEInit/D");
  simTree.Branch("endpointKE_u", &endpointKE_u, "endpointKE_u/D");
  simTree.Branch("endpointKE_p", &endpointKE_p, "endpointKE_p/D");
  simTree.Branch("endpointPitchAngle", &endpointPitchAngle,
                 "endpointPitchAngle/D");
  simTree.Branch("perturberStartPos", perturberStartPos,
                 "perturberStartPos[3]/D");
  simTree.Branch("perturberStartVel", perturberStartVel,
                 "perturberStartVel[3]/D");
  simTree.Branch("perturberKEInit", &perturberKEInit, "perturberKEInit/D");
  simTree.Branch("perturberPitchAngle", &perturberPitchAngle,
                 "perturberPitchAngle/D");

  while (simCounter < nSims) {
    std::cout << "Simulation " << simCounter << std::endl;
    // Start a clock
    std::clock_t start{std::clock()};

    TString trackFileNoPerturbation{outputDir + "/track_" + make_uuid() +
                                    ".root"};
    TVector3 posGenEndpoint(0, 0, 0);
    TVector3 velGenEndpoint(0, 0, 0);
    double EFinalUnperturbed{0};
    double tFinalUnperturbed{0};

    // Let's generate an endpoint electron
    bool GoodSim{false};
    while (!GoodSim) {
      std::random_device rd;
      std::mt19937 mt(rd());
      std::uniform_real_distribution<double> uniDist(0, 1);
      // Generate a random position and direction in the trap
      const double genLength{0.5 * trapLength};
      double zGen{uniDist(mt) * genLength - genLength / 2};
      double rhoGen{sqrt(uniDist(mt)) * trapCoilRadius};
      double phiPosGen{uniDist(mt) * 2 * M_PI};
      posGenEndpoint.SetX(rhoGen * cos(phiPosGen));
      posGenEndpoint.SetY(rhoGen * sin(phiPosGen));
      posGenEndpoint.SetZ(zGen);
      velGenEndpoint = GenerateRandomDirection(mt) * endpointSpeed;

      // Get the pitch angle
      double pitchAngle{abs(atan(velGenEndpoint.Perp() / velGenEndpoint.Z()))};
      if (pitchAngle > thetaMin) {
        std::cout << "Start position = " << posGenEndpoint.X() << " "
                  << posGenEndpoint.Y() << " " << posGenEndpoint.Z()
                  << std::endl;
        std::cout << "Pitch angle = " << pitchAngle * 180 / M_PI
                  << " degrees\n";

        endpointKEInit = endpointEnergy;
        endpointPitchAngle = pitchAngle;
        endpointStartPos[0] = posGenEndpoint.X();
        endpointStartPos[1] = posGenEndpoint.Y();
        endpointStartPos[2] = posGenEndpoint.Z();
        endpointStartVel[0] = velGenEndpoint.X();
        endpointStartVel[1] = velGenEndpoint.Y();
        endpointStartVel[2] = velGenEndpoint.Z();

        // Run the simulation
        std::cout << "Writing track file to " << trackFileNoPerturbation
                  << std::endl;
        ElectronTrajectoryGen traj(trackFileNoPerturbation, bathtubField,
                                   posGenEndpoint, velGenEndpoint, dt, tMax,
                                   0.0, 2 * R_E / (3 * TMath::C()));

        // Open the file and get the last entry of the TTree
        TFile tFileNoPerturbation(trackFileNoPerturbation, "READ");
        TTree* trajTree = (TTree*)tFileNoPerturbation.Get("tree");
        double vx, vy, vz;
        double t;
        trajTree->SetBranchAddress("time", &t);
        trajTree->SetBranchAddress("xVel", &vx);
        trajTree->SetBranchAddress("yVel", &vy);
        trajTree->SetBranchAddress("zVel", &vz);
        trajTree->GetEntry(trajTree->GetEntries() - 1);
        double beta{sqrt(vx * vx + vy * vy + vz * vz) / TMath::C()};
        EFinalUnperturbed = (1 / sqrt(1 - beta * beta)) * ME_EV;
        endpointKE_u = EFinalUnperturbed - ME_EV;

        tFinalUnperturbed = t;
        GoodSim = true;
        tFileNoPerturbation.Close();
        remove(trackFileNoPerturbation);
      }
    }

    // Now we need to generate a second electron to perturb the first
    bool GoodSim2{false};
    TVector3 posGenPerturber(0, 0, 0);
    TVector3 velGenPerturber(0, 0, 0);
    while (!GoodSim2) {
      std::random_device rd;
      std::mt19937 mt(rd());
      std::uniform_real_distribution<double> uniDist(0, 1);
      // Start with the endpoint electron
      // Generate a random position in the trap
      const double genLength{0.5 * trapLength};
      double zGen{uniDist(mt) * genLength - genLength / 2};
      double rhoGen{sqrt(uniDist(mt)) * trapCoilRadius};
      double phiPosGen{uniDist(mt) * 2 * M_PI};
      posGenPerturber.SetX(rhoGen * cos(phiPosGen));
      posGenPerturber.SetY(rhoGen * sin(phiPosGen));
      posGenPerturber.SetZ(zGen);

      // Depending on selected option, draw a velocity from the spectrum or take
      // the most average value of the spectrum
      const double engAvg{5.7e3};
      const double speedAvg{GetSpeedFromKE(engAvg, ME)};
      velGenPerturber = (drawFromSpectrum)
                            ? GenerateDecayVelocity(p, mt)
                            : GenerateRandomDirection(mt) * speedAvg;

      // Get the pitch angle
      double pitchAngle{
          abs(atan(velGenPerturber.Perp() / velGenPerturber.Z()))};

      if (pitchAngle > thetaMin) {
        std::cout << "Start position = " << posGenPerturber.X() << " "
                  << posGenPerturber.Y() << " " << posGenPerturber.Z()
                  << std::endl;
        std::cout << "Pitch angle = " << pitchAngle * 180 / M_PI
                  << " degrees\n";

        perturberStartPos[0] = posGenPerturber.X();
        perturberStartPos[1] = posGenPerturber.Y();
        perturberStartPos[2] = posGenPerturber.Z();
        perturberStartVel[0] = velGenPerturber.X();
        perturberStartVel[1] = velGenPerturber.Y();
        perturberStartVel[2] = velGenPerturber.Z();
        perturberPitchAngle = pitchAngle;
        perturberKEInit =
            ((1 / sqrt(1 - pow(velGenPerturber.Mag() / TMath::C(), 2))) - 1) *
            ME_EV;

        GoodSim2 = true;
      }
    }

    // Now let's use the multiparticle solver to simulate the two electrons
    Particle original(posGenEndpoint, velGenEndpoint);
    Particle perturber(posGenPerturber, velGenPerturber);
    std::vector<Particle> twoParticleListPerturber{original, perturber};
    ManyParticleSolver twoElectronSolverPerturber(bathtubField,
                                                  twoParticleListPerturber);
    double currentTime{0};
    int nSteps{0};
    while (currentTime <= tFinalUnperturbed) {
      twoElectronSolverPerturber.AdvanceParticles(dt);
      nSteps++;
      currentTime = double(nSteps) * dt;
    }

    // Having finished running, get the final energy of the endpoint electron
    TVector3 finalVel{twoElectronSolverPerturber.GetParticle(0).GetVelocity()};
    double beta{finalVel.Mag() / TMath::C()};
    double EFinalPerturbed{(1 / sqrt(1 - beta * beta)) * ME_EV};
    endpointKE_p = EFinalPerturbed - ME_EV;
    simTree.Fill();

    // Print the final energies
    std::cout
        << "Final energy of electrons (unperturbed, perturbed, difference) = "
        << EFinalUnperturbed << " eV,\t" << EFinalPerturbed << " eV,\t"
        << (EFinalPerturbed - EFinalUnperturbed) * 1e3 << " meV\n";

    simCounter++;
    std::clock_t end{std::clock()};
    // Print the time taken in seconds
    std::cout << "Time taken: " << double(end - start) / CLOCKS_PER_SEC
              << " seconds\n";
  }

  outFile.cd();
  simTree.Write();
  outFile.Close();
  return 0;
}