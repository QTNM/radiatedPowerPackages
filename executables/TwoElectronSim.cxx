/*
  TwoElectronSim.cxx

  Test program to look at the motion of two electrons in a magnetic trap.
*/

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

TVector3 GenerateRandomDirection(std::mt19937& mt) {
  std::uniform_real_distribution<double> uniDist(0, 1);
  double phiGen{uniDist(mt) * 2 * M_PI};
  double thetaGen{uniDist(mt) * M_PI};
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
  TString outputFile{argv[1]};
  TFile outFile(outputFile, "RECREATE");

  // We want to look at the effect on an endpoint electron
  const double endpointEnergy{18.6e3};                             // eV
  const double endpointSpeed{GetSpeedFromKE(endpointEnergy, ME)};  // m/s

  // Simulation timings setup
  const double dt{1e-12};    // s
  const double tMax{10e-6};  // s

  // Magnetic field setup
  const double B0{1.0};                                                  // T
  const double trapDepth{4e-3};                                          // T
  const double trapCoilRadius{0.03};                                     // m
  const double trapCoilCurrent{2.0 * trapDepth * trapCoilRadius / MU0};  // A
  const double trapLength{0.2};                                          // m
  auto bathtubField{new BathtubField(trapCoilRadius, trapCoilCurrent,
                                     -trapLength / 2, trapLength / 2,
                                     TVector3(0, 0, B0))};

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
  while (!GoodSim) {
    std::random_device rd;
    std::mt19937 mt(rd());
    std::uniform_real_distribution<double> uniDist(0, 1);
    // Start with the endpoint electron
    // Generate a random position in the trap
    double zGen{uniDist(mt) * trapLength - trapLength / 2};
    double rhoGen{uniDist(mt) * trapCoilRadius};
    double phiPosGen{uniDist(mt) * 2 * M_PI};
    TVector3 posGen(rhoGen * cos(phiPosGen), rhoGen * sin(phiPosGen), zGen);
    // Generate a random direction and use the speed
    TVector3 velGen{GenerateRandomDirection(mt) * endpointSpeed};

    // Now generate the second electron
    TVector3 velGen2{GenerateDecayVelocity(p, mt)};
    // Generate a second random position
    double zGen2{uniDist(mt) * trapLength - trapLength / 2};
    double rhoGen2{uniDist(mt) * trapCoilRadius};
    double phiPosGen2{uniDist(mt) * 2 * M_PI};
    TVector3 posGen2(rhoGen2 * cos(phiPosGen2), rhoGen2 * sin(phiPosGen2),
                     zGen2);
  }

  outFile.Close();
  return 0;
}