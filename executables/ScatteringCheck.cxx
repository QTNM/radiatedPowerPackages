/*
  ScatteringCheck.cxx
*/

#include <getopt.h>
#include <unistd.h>

#include <cmath>
#include <iostream>
#include <memory>
#include <random>
#include <vector>

#include "BasicFunctions/BasicFunctions.h"
#include "BasicFunctions/Constants.h"
#include "BasicFunctions/TritiumSpectrum.h"
#include "ElectronDynamics/BorisSolver.h"
#include "ElectronDynamics/QTNMFields.h"
#include "Scattering/ElasticScatter.h"
#include "Scattering/InelasticScatter.h"
#include "TBranch.h"
#include "TFile.h"
#include "TH1.h"
#include "TString.h"
#include "TTree.h"
#include "TVector3.h"

using namespace rad;
using std::cout;
using std::endl;
using std::vector;

typedef std::piecewise_linear_distribution<double> pld;

pld GenSpectrum(int nPnts) {
  vector<double> eVec{};
  vector<double> rateVec{};
  const double endpoint{18575};

  for (int i{0}; i < nPnts; i++) {
    double E{0.1 + endpoint * double(i) / double(nPnts - 1)};
    eVec.push_back(E);
    // Assume massless neutrinos
    rateVec.push_back(TritiumDecayRate(E, 0, 0, 0, endpoint));
  }

  return pld(eVec.begin(), eVec.end(), rateVec.begin());
}

/// @brief Generates a random velocity from a inputted ke distribution
/// @param p PDF of kinetic energies
/// @param mt Random number generator
/// @return A velocity vector in units of m s^-1
TVector3 GenDecayVelocity(pld &p, std::mt19937 &mt) {
  double ke{p(mt)};  // Kinetic energy in eV

  // Calculate electron speed
  double v{GetSpeedFromKE(ke, ME)};

  // Get a random direction
  std::uniform_real_distribution<double> uniDist(0, 1);
  double phiGen{uniDist(mt) * TMath::TwoPi()};
  double thetaGen{uniDist(mt) * TMath::Pi()};

  return TVector3(v * cos(phiGen) * sin(thetaGen),
                  v * sin(phiGen) * sin(thetaGen), v * cos(thetaGen));
}

TVector3 RotateToCoords(TVector3 v, TVector3 newX, TVector3 newY,
                        TVector3 newZ) {
  // We are just transforming from the ROOT frame so our old axis coordinates
  // are just the unit vectors
  TVector3 oldX(1, 0, 0);
  TVector3 oldY(0, 1, 0);
  TVector3 oldZ(0, 0, 1);
  double pXPrime{newX.X() * v.X() + newY.X() * v.Y() + newZ.X() * v.Z()};
  double pYPrime{newX.Y() * v.X() + newY.Y() * v.Y() + newZ.Y() * v.Z()};
  double pZPrime{newX.Z() * v.X() + newY.Z() * v.Y() + newZ.Z() * v.Z()};
  TVector3 newVector(pXPrime, pYPrime, pZPrime);
  newVector = newVector.Unit();
  return newVector;
}

/// @brief Reset an array to nonsense
/// @param arr The array in question
/// @param arrSize Size of array
void ResetArray(unsigned int *arr, unsigned int arrSize) {
  for (unsigned int i{0}; i < arrSize; i++) arr[i] = 0;
}

/// @brief Reset an array to nonsense
/// @param arr The array in question
/// @param arrSize Size of array
void ResetArray(double *arr, unsigned int arrSize) {
  for (unsigned int i{0}; i < arrSize; i++) arr[i] = -1;
}

int main(int argc, char *argv[]) {
  int opt;

  std::string outputFile{" "};
  unsigned int nElectrons{100};

  while ((opt = getopt(argc, argv, ":o:n:")) != -1) {
    switch (opt) {
      case 'o':
        outputFile = optarg;
        std::cout << "Output file is " << outputFile << std::endl;
        break;

      case 'n':
        nElectrons = std::stoi(optarg);
        break;

      case ':':
        std::cout << "Option needs a value\n";
        break;

      case '?':
        std::cout << "Unknown option: " << optopt << std::endl;
        break;
    }
  }

  // Check input
  if (outputFile == " ") {
    std::cout << "Must specify output file with -o" << std::endl;
    exit(1);
  }

  std::cout << "Simulating " << nElectrons << " electrons\n";

  TFile fout(outputFile.data(), "recreate");

  std::random_device rd;
  std::mt19937 gen(rd());

  // First generate the spectrum distribution
  const int nSpecPnts{2000};
  pld spec{GenSpectrum(nSpecPnts)};

  // Create TTree with the same information in (easier to combine)
  TTree outTree("scatTree", "Scattering tree");
  const unsigned int nMaxScatters{30};

  double eInit{0};
  double totalTime{0};
  unsigned int nScatters{0};
  double scatterLen[nMaxScatters];
  double scatterTime[nMaxScatters];
  double scatterAng[nMaxScatters];
  unsigned int scatterEl[nMaxScatters];
  unsigned int scatterInel[nMaxScatters];
  double scatterELoss[nMaxScatters];

  outTree.Branch("eInit", &eInit, "eInit/D");
  outTree.Branch("nScatters", &nScatters, "nScatters/i");
  outTree.Branch("lifetime", &totalTime, "lifetime/D");
  outTree.Branch("scatterLen", scatterLen, "scatterLen[nScatters]/D");
  outTree.Branch("scatterTime", scatterTime, "scatterTime[nScatters]/D");
  outTree.Branch("scatterAng", scatterAng, "scatterAng[nScatters]/D");
  outTree.Branch("scatterEl", scatterEl, "scatterEl[nScatters]/i");
  outTree.Branch("scatterInel", scatterInel, "scatterInel[nScatters]/i");
  outTree.Branch("scatterELoss", scatterELoss, "scatterELoss[nScatters]/D");

  // Define a magnetic trap. Start with a harmonic one
  const double coilRadius{0.03};                               // m
  const double trapDepth{4e-3};                                // T
  const double coilCurrent{2 * trapDepth * coilRadius / MU0};  // T
  const double bkg{0.72};                                      // Tesla
  auto field = new HarmonicField(coilRadius, coilCurrent, bkg);
  const double tritiumDensity{2e18};                           // atoms m^-3
  const double zMax{0.1};                                      // metres

  for (unsigned int i{0}; i < nElectrons; i++) {
    TVector3 pos(0, 0, 0);  // metres
    TVector3 vel{GenDecayVelocity(spec, gen)};
    double speed{vel.Mag()};
    // Recover the electron kinetic energy
    double gamma{1 / sqrt(1 - pow(vel.Mag() / TMath::C(), 2))};
    double ke{(gamma - 1) * ME_EV};

    cout << "Electron " << i + 1 << ": E_i = " << ke / 1e3 << " keV\n";

    bool isTrapped{true};
    const double keCutoff{1};       // eV
    const double timeCutoff{1e-3};  // s
    const double tau{2 * R_E / (3 * TMath::C())};

    // Variables used in the tree
    totalTime = 0;  // s
    eInit = ke;
    nScatters = 0;
    ResetArray(scatterLen, nMaxScatters);
    ResetArray(scatterTime, nMaxScatters);
    ResetArray(scatterAng, nMaxScatters);
    ResetArray(scatterEl, nMaxScatters);
    ResetArray(scatterInel, nMaxScatters);
    ResetArray(scatterELoss, nMaxScatters);

    double scatAngle{0};
    while (isTrapped && ke > 1 && totalTime < timeCutoff &&
           nScatters < nMaxScatters) {
      // How far does this electron travel before its first scatter?
      ElasticScatter scatEl(ke);
      double elXSec{scatEl.GetTotalXSec()};
      InelasticScatter scatInel(ke);
      double inelXSec{scatInel.GetTotalXSec()};
      double totalXSec{elXSec + inelXSec};

      // Calculate the mean free path
      const double lambdaStep{1 / (tritiumDensity * totalXSec)};
      std::exponential_distribution<double> pathDistStep(1 / lambdaStep);
      const double pathLenStep{pathDistStep(gen)};
      const double pathTimeStep{pathLenStep / vel.Mag()};

      scatterLen[nScatters] = pathLenStep;
      scatterTime[nScatters] = pathTimeStep;

      std::cout << "Scattering after " << pathTimeStep * 1e6 << " us\n";
      // Propagate the particle up to its next scatter
      // Make sure to check that it has not exited the trap
      double simStepTime{5e-12};  // seconds
      BorisSolver solver(field, -TMath::Qe(), ME, tau);

      double nTimeSteps{std::round(pathTimeStep / simStepTime)};
      for (int iStep{0}; iStep < nTimeSteps; iStep++) {
        // Advance time
        totalTime += simStepTime;
        if (std::fmod(totalTime, 5e-6) < simStepTime) {
          cout << "Simulated " << totalTime * 1e6 << " us\n";
        }

        auto outputStep = solver.advance_step(simStepTime, pos, vel);
        // Update position and velocity vectors
        pos = std::get<0>(outputStep);
        vel = std::get<1>(outputStep);

        // Check if we are out of bounds
        if (abs(pos.Z()) > zMax) {
          isTrapped = false;
          break;
        }
      }  // Loop over propagation steps

      if (isTrapped) {
        // Now figure out kinematics of the scatter
        // Calculate the current kinetic energy
        gamma = 1 / sqrt(1 - pow(vel.Mag() / TMath::C(), 2));
        ke = (gamma - 1) * ME_EV;

        // Recalculate the cross sections based on the cross-sections
        ElasticScatter scatEl2(ke);
        InelasticScatter scatInel2(ke);
        elXSec = scatEl2.GetTotalXSec();
        inelXSec = scatInel2.GetTotalXSec();
        totalXSec = elXSec + inelXSec;
        cout << "Prob (elastic, inelastic) = (" << elXSec / totalXSec << ", "
             << inelXSec / totalXSec << ")\n";

        // First figure out if this is an elastic or inelastic scatter
        std::uniform_real_distribution<double> uni(0.0, 1.0);
        double eLoss{0};
        if (uni(gen) < elXSec / totalXSec) {
          // We have an elastic scatter
          cout << "Elastic scatter\n";
          // No energy loss so just get the scattering angle
          scatAngle = scatEl2.GetRandomScatteringAngle();

          scatterEl[nScatters] = 1;
          scatterInel[nScatters] = 0;
        } else {
          // We have an inelastic scatter
          cout << "Inelastic scatter\n";
          // Get the energy of the secondary
          double wSample{scatInel2.GetRandomW()};
          double theta2Sample{scatInel2.GetRandomTheta(wSample)};
          // Now calculate the energy and the scattering angle of the primary
          scatAngle = scatInel2.GetPrimaryScatteredAngle(wSample, theta2Sample);

          eLoss = ke - scatInel2.GetPrimaryScatteredE(wSample, theta2Sample);
          ke = scatInel2.GetPrimaryScatteredE(wSample, theta2Sample);

          scatterEl[nScatters] = 0;
          scatterInel[nScatters] = 1;
        }

        speed = GetSpeedFromKE(ke, ME);

        scatterAng[nScatters] = scatAngle;
        scatterELoss[nScatters] = eLoss;

        cout << "Scattering angle = " << scatAngle * 180 / TMath::Pi()
             << " degrees\t New KE = " << ke / 1e3
             << " keV\t Energy loss = " << eLoss << " eV\n";

        // Now we have to calculate the new direction of the electron
        // We have the polar angle already from sampling
        // Azimuthal angle distributed uniformly
        double phiSampled{uni(gen) * TMath::TwoPi()};
        TVector3 originalDir{vel.Unit()};
        // Define the new direction in the frame of the particle
        TVector3 oldDir(0, 0, 1);
        TVector3 newDir(sin(scatAngle) * cos(phiSampled),
                        sin(scatAngle) * sin(phiSampled), cos(scatAngle));
        TVector3 ax2(-originalDir.Y() / originalDir.X(), 1, 0);
        ax2 = ax2.Unit();
        TVector3 ax3{(ax2.Cross(originalDir)).Unit()};
        // Now you can get the third axis we want to rotate to
        TVector3 newDirection{RotateToCoords(newDir, ax2, ax3, originalDir)};

        // Set the new velocity
        vel = speed * newDirection;

        // Move onto the next scatter
        nScatters++;
      }
    }

    // We have now completed the motion
    cout << "Escaped after " << totalTime * 1e6 << " us and " << nScatters
         << " scatters\n\n";

    outTree.Fill();
  }
  delete field;

  fout.cd();
  outTree.Write();

  fout.Close();
  return 0;
}