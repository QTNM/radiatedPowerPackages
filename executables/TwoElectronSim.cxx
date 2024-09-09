/*
  TwoElectronSim.cxx

  Test program to look at the motion of two electrons in a magnetic trap.
*/

#include <getopt.h>
#include <unistd.h>

#include <iostream>
#include <memory>
#include <random>
#include <vector>

#include "BasicFunctions/BasicFunctions.h"
#include "BasicFunctions/Constants.h"
#include "BasicFunctions/TritiumSpectrum.h"
#include "ElectronDynamics/BorisSolver.h"
#include "ElectronDynamics/QTNMFields.h"
#include "ElectronDynamics/TrajectoryGen.h"
#include "TFile.h"
#include "TGraph.h"
#include "TString.h"
#include "TTree.h"
#include "TVector3.h"

using namespace rad;

typedef std::piecewise_linear_distribution<double> pld;

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
  std::string outputFile;
  bool drawFromSpectrum{false};
  int opt{};
  while ((opt = getopt(argc, argv, ":o:r")) != -1) {
    switch (opt) {
      case 'o':
        std::cout << "Output file: " << optarg << std::endl;
        outputFile = optarg;
        break;
      case 'r':
        std::cout << "Drawing from spectrum" << std::endl;
        drawFromSpectrum = true;
        break;
      case ':':
        std::cerr << "Option -" << optopt << " requires an argument.\n";
        return 1;
      case '?':
        std::cerr << "Unknown option -" << optopt << std::endl;
        return 1;
    }
  }

  TFile outFile(outputFile.c_str(), "RECREATE");

  // We want to look at the effect on an endpoint electron
  const double endpointEnergy{18.575e3};                           // eV
  const double endpointSpeed{GetSpeedFromKE(endpointEnergy, ME)};  // m/s

  // Simulation timings setup
  const double dt{1e-12};    // s
  const double tMax{10e-6};  // s

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

  bool GoodSim{false};
  pld p{GenerateSpectrum(1000)};

  TString trackFile{"track.root"};
  while (!GoodSim) {
    std::random_device rd;
    std::mt19937 mt(rd());
    std::uniform_real_distribution<double> uniDist(0, 1);
    // Start with the endpoint electron
    // Generate a random position in the trap
    const double genLength{0.5 * trapLength};
    double zGen{uniDist(mt) * genLength - genLength / 2};
    double rhoGen{uniDist(mt) * trapCoilRadius};
    double phiPosGen{uniDist(mt) * 2 * M_PI};
    TVector3 posGen(rhoGen * cos(phiPosGen), rhoGen * sin(phiPosGen), zGen);
    // Depending on selected option, draw a velocity from the spectrum of take
    // the most probable value of the spectrum
    const double engAvg{5.7e3};
    const double speedAvg{GetSpeedFromKE(engAvg, ME)};
    TVector3 velGen = (drawFromSpectrum)
                          ? GenerateDecayVelocity(p, mt)
                          : GenerateRandomDirection(mt) * speedAvg;

    // Get the pitch angle
    double pitchAngle{abs(atan(velGen.Perp() / velGen.Z()))};

    if (pitchAngle > thetaMin) {
      std::cout << "Start position = " << posGen.X() << " " << posGen.Y() << " "
                << posGen.Z() << std::endl;
      std::cout << "Pitch angle = " << pitchAngle * 180 / M_PI << " degrees\n";

      // Run the simulation
      std::cout << "Writing track file to " << trackFile << std::endl;
      ElectronTrajectoryGen traj(trackFile, bathtubField, posGen, velGen, dt,
                                 tMax, 0.0, 2 * R_E / (3 * TMath::C()));
      GoodSim = true;
    }
  }

  // Now let's generate an endpoint electron
  bool GoodSim2{false};
  while (!GoodSim2) {
    std::random_device rd;
    std::mt19937 mt(rd());
    std::uniform_real_distribution<double> uniDist(0, 1);
    // Generate a random position and direction in the trap
    const double genLength{0.5 * trapLength};
    double zGen{uniDist(mt) * genLength - genLength / 2};
    double rhoGen{uniDist(mt) * trapCoilRadius};
    double phiPosGen{uniDist(mt) * 2 * M_PI};
    TVector3 posGen(rhoGen * cos(phiPosGen), rhoGen * sin(phiPosGen), zGen);
    TVector3 velGen{GenerateRandomDirection(mt) * endpointSpeed};

    // Get the pitch angle
    double pitchAngle{abs(atan(velGen.Perp() / velGen.Z()))};
    if (pitchAngle > thetaMin) {
      std::cout << "Start position = " << posGen.X() << " " << posGen.Y() << " "
                << posGen.Z() << std::endl;
      std::cout << "Pitch angle = " << pitchAngle * 180 / M_PI << " degrees\n";
      GoodSim2 = true;
    }
  }

  // Delete the first track file
  remove(trackFile.Data());

  outFile.Close();
  return 0;
}