/*
  SecondaryElectronProduction.cxx

  Executable checking the trapping time of primary electrons along with types of
  secondary electrons produced.
*/

#include <getopt.h>
#include <unistd.h>

#include <boost/lexical_cast.hpp>
#include <boost/uuid/uuid.hpp>
#include <boost/uuid/uuid_generators.hpp>
#include <boost/uuid/uuid_io.hpp>
#include <cmath>
#include <ctime>
#include <iostream>
#include <random>
#include <string>

#include "BasicFunctions/BasicFunctions.h"
#include "BasicFunctions/Constants.h"
#include "ElectronDynamics/BorisSolver.h"
#include "ElectronDynamics/QTNMFields.h"
#include "Scattering/ElasticScatter.h"
#include "Scattering/InelasticScatter.h"
#include "TFile.h"
#include "TString.h"
#include "TTree.h"
#include "TVector3.h"

using namespace rad;
using std::cout;
using std::endl;

/// @brief Function for generating a random file string
/// @return A unique string which can be used for file names
std::string make_uuid() {
  return boost::lexical_cast<std::string>((boost::uuids::random_generator())());
}

// Main program
int main(int argc, char* argv[]) {
  int opt{};
  std::string outputDirName{" "};
  unsigned int nSims{100};
  double tritiumDensity{1e18};  // atoms/m^3

  while ((opt = getopt(argc, argv, "o:n:d:")) != -1) {
    switch (opt) {
      case 'o':
        outputDirName = optarg;
        break;
      case 'n':
        nSims = std::stoi(optarg);
        break;
      case 'd':
        tritiumDensity = std::stod(optarg);
        break;
      case ':':
        std::cerr << "Option -" << optopt << " requires an argument."
                  << std::endl;
        return 1;
      case '?':
        std::cerr << "Unrecognised option: -" << optopt << std::endl;
        return 1;
      default:
        std::cerr << "Usage: " << argv[0] << " [-o output directory] [-n nSims]"
                  << std::endl;
        return 1;
    }
  }

  cout << "Output directory: " << outputDirName << endl;
  cout << "Attempting to generate " << nSims << " simulations." << endl;
  cout << "Tritium density: " << tritiumDensity << " atoms/m^3" << endl;

  // Create the output file
  TFile outFile(TString::Format("%s/secElec_%s.root", outputDirName.c_str(),
                                make_uuid().c_str()),
                "RECREATE");
  TTree outTree("secElecTree", "Secondary electron production");
  outTree.SetDirectory(&outFile);
  const size_t nMaxScatters = 100;
  int nScatters;
  double exitTime;
  outTree.Branch("nScatters", &nScatters, "nScatters/I");
  outTree.Branch("exitTime", &exitTime, "exitTime/D");

  double startPos[3];
  double startVel[3];
  double startKE;
  outTree.Branch("startPos", startPos, "startPos[3]/D");
  outTree.Branch("startVel", startVel, "startVel[3]/D");
  outTree.Branch("startKE", &startKE, "startKE/D");

  double scatterTime[nMaxScatters];
  double incidentPosX[nMaxScatters];
  double incidentPosY[nMaxScatters];
  double incidentPosZ[nMaxScatters];
  double incidentVelX[nMaxScatters];
  double incidentVelY[nMaxScatters];
  double incidentVelZ[nMaxScatters];
  outTree.Branch("scatterTime", scatterTime, "scatterTime[nScatters]/D");
  outTree.Branch("incidentPosX", incidentPosX, "incidentPosX[nScatters]/D");
  outTree.Branch("incidentPosY", incidentPosY, "incidentPosY[nScatters]/D");
  outTree.Branch("incidentPosZ", incidentPosZ, "incidentPosZ[nScatters]/D");
  outTree.Branch("incidentVelX", incidentVelX, "incidentVelX[nScatters]/D");
  outTree.Branch("incidentVelY", incidentVelY, "incidentVelY[nScatters]/D");
  outTree.Branch("incidentVelZ", incidentVelZ, "incidentVelZ[nScatters]/D");

  double scatterAngle[nMaxScatters];
  double incidentKE[nMaxScatters];
  double secElecKE[nMaxScatters];
  outTree.Branch("scatterAngle", scatterAngle, "scatterAngle[nScatters]/D");
  outTree.Branch("incidentKE", incidentKE, "incidentKE[nScatters]/D");
  outTree.Branch("secElecKE", secElecKE, "secElecKE[nScatters]/D");

  // Define a big bathtub trap with the electrons
  const double trapRadius{0.3};
  const double trapLength{3.5};
  const double trappingFraction{0.1};
  const double deltaTheta{asin(trappingFraction)};
  const double BMin{1.0};                                          // Tesla
  const double BMax{BMin / pow(cos(deltaTheta), 2)};               // Tesla
  const double coilCurrent{2 * (BMax - BMin) * trapRadius / MU0};  // Amps
  auto field{new BathtubField(trapRadius, coilCurrent, -trapLength / 2,
                              trapLength / 2, TVector3(0, 0, BMin))};

  // Set up random number stuff
  std::random_device rd{};
  std::mt19937 gen(rd());

  // Calc the step size we want here
  const double EMax{18.6e3};  // eV
  const double EMin{100.0};   // eV
  const double maxCycFreq{CalcCyclotronFreq(EMax, BMin)};
  const double deltaT{1 / (10 * maxCycFreq)};  // seconds
  cout << "Time step: " << deltaT << " seconds" << endl;

  const double tau{2 * R_E / (3 * TMath::C())};
  for (unsigned int iEv{0}; iEv < nSims; iEv++) {
    cout << "Simulating event " << iEv + 1 << endl;

    // Reset the branch variables
    nScatters = 0;
    exitTime = 0;
    for (size_t i{0}; i < nMaxScatters; i++) {
      scatterTime[i] = -1;
      incidentPosX[i] = 0;
      incidentPosY[i] = 0;
      incidentPosZ[i] = 0;
      incidentVelX[i] = 0;
      incidentVelY[i] = 0;
      incidentVelZ[i] = 0;
      scatterAngle[i] = -1;
      incidentKE[i] = -1;
      secElecKE[i] = -1;
    }

    // Generate a random position in the cylinder
    std::uniform_real_distribution<double> uni1(0, 1);
    double zGen{trapLength * uni1(gen) - trapLength / 2};
    double thetaPosGen{uni1(gen) * 2 * M_PI};
    double rGen{trapRadius * sqrt(uni1(gen))};

    TVector3 genPos(rGen * cos(thetaPosGen), rGen * sin(thetaPosGen), zGen);

    // Generate a random energy uniform across the spectrum
    double EGen{EMin + (EMax - EMin) * uni1(gen)};
    const double vGen{GetSpeedFromKE(EGen, ME)};
    // Generate isotropic velocity
    const double phiVelGen{uni1(gen) * 2 * M_PI};
    const double thetaVelGen{acos(2 * uni1(gen) - 1)};
    TVector3 genVel(vGen * sin(thetaVelGen) * cos(phiVelGen),
                    vGen * sin(thetaVelGen) * sin(phiVelGen),
                    vGen * cos(thetaVelGen));
    startKE = EGen;
    startPos[0] = genPos.X();
    startPos[1] = genPos.Y();
    startPos[2] = genPos.Z();
    startVel[0] = genVel.X();
    startVel[1] = genVel.Y();
    startVel[2] = genVel.Z();

    cout << "Initial energy: " << EGen << " eV" << endl;
    cout << "Initial position: (" << genPos.X() << ", " << genPos.Y() << ", "
         << genPos.Z() << ") m" << endl;
    cout << "Initial velocity: (" << genVel.X() << ", " << genVel.Y() << ", "
         << genVel.Z() << ") m/s" << endl;

    double tSim{0};
    const double maxTSim{1e-3};    // seconds
    const double printTime{10e-6};  // seconds
    // Set up the Boris solver
    BorisSolver solver(field, -TMath::Qe(), ME, tau);

    while (tSim < maxTSim && abs(genPos.Z()) < trapLength / 2) {
      // Set up print out for the time
      if (fmod(tSim, printTime) < deltaT) {
        cout << "Time: " << tSim * 1e6 << " us\tEnergy = " << EGen << " eV"
             << endl;
      }

      // Set up the cross sections
      ElasticScatter elasticScatter(EGen);
      InelasticScatter inelasticScatter(EGen);
      const double elasticXSec{elasticScatter.GetTotalXSec()};
      const double inelasticXSec{inelasticScatter.GetTotalXSec()};
      const double totalXSec{elasticXSec + inelasticXSec};
      const double meanFreeTime{1 /
                                (totalXSec * tritiumDensity * genVel.Mag())};
      // Draw thew scattering time from an exponential distribution
      std::exponential_distribution<double> expDist(1 / meanFreeTime);
      const double tScatter{expDist(gen)};
      // Check if the scattering time is less than the time step
      if (tScatter < deltaT) {
        cout << "Scattering event at " << (tSim + tScatter) * 1e6 << " us"
             << endl;
        tSim += tScatter;
        auto [pos, vel] = solver.advance_step(tScatter, genPos, genVel);
        genPos = pos;
        genVel = vel;
        EGen = (sqrt(1 / (1 - pow(genVel.Mag() / TMath::C(), 2))) - 1) * ME_EV;

        scatterTime[nScatters] = tSim;
        incidentPosX[nScatters] = genPos.X();
        incidentPosY[nScatters] = genPos.Y();
        incidentPosZ[nScatters] = genPos.Z();
        incidentVelX[nScatters] = genVel.X();
        incidentVelY[nScatters] = genVel.Y();
        incidentVelZ[nScatters] = genVel.Z();
        incidentKE[nScatters] = EGen;

        // Check if this is an inelastic or elastic scattering event
        if (uni1(gen) < elasticXSec / totalXSec) {
          // Elastic scattering
          const double theta{elasticScatter.GetRandomScatteringAngle()};
          cout << "Elastic scattering event" << endl;
          cout << "Scattering angle: " << theta * 180 / M_PI << " degrees"
               << endl;

          scatterAngle[nScatters] = theta;

          // Update the velocity
          genVel = elasticScatter.GetScatteredVector(genVel, EGen, theta);
        } else {
          // Inelastic scattering
          // Get the energy of the outgoing electron
          const double W{inelasticScatter.GetRandomW()};
          const double theta{inelasticScatter.GetRandomTheta(W)};

          scatterAngle[nScatters] = theta;
          secElecKE[nScatters] = EGen - W - RYDBERG_EV;

          cout << "Inelastic scattering event" << endl;
          cout << "Scattering angle: " << theta * 180 / M_PI << " degrees"
               << endl;
          cout << "Outgoing energy: " << W << " eV" << endl;
          // Update the energy and velocity
          genVel = inelasticScatter.GetScatteredVector(genVel, W, theta);
          EGen = W;
        }
        nScatters++;
      } else {
        auto [pos, vel] = solver.advance_step(deltaT, genPos, genVel);
        tSim += deltaT;
        genPos = pos;
        genVel = vel;
        EGen = (sqrt(1 / (1 - pow(genVel.Mag() / TMath::C(), 2))) - 1) * ME_EV;
      }
    }

    cout << "Simulation ended at " << tSim * 1e6 << " us\n" << endl;
    exitTime = tSim;
    outTree.Fill();
  }

  outFile.cd();
  outTree.Write();

  outFile.Close();
  delete field;
  return 0;
}
