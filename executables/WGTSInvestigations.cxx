/*
  WGTSInvestigations.cxx

  Look at some of the numbers from Nastoyashchii, Titov, Morozov, Gluck & Otten
  (2005) and see if they can be verified.
*/

#include <getopt.h>
#include <unistd.h>

#include <cmath>
#include <iostream>
#include <random>
#include <tuple>
#include <vector>

#include "BasicFunctions/Constants.h"
#include "BasicFunctions/TritiumSpectrum.h"
#include "ElectronDynamics/BorisSolver.h"
#include "ElectronDynamics/QTNMFields.h"
#include "Scattering/ElasticScatter.h"
#include "Scattering/InelasticScatter.h"
#include "TFile.h"
#include "TMath.h"
#include "TTree.h"
#include "TVector3.h"

using namespace rad;

typedef std::piecewise_linear_distribution<double> pld;

pld GenerateSpectrum(unsigned int nPoints) {
  std::vector<double> energyVec{};
  std::vector<double> decayRateVec{};
  const double endpointEnergy{18.575e3};

  for (unsigned int i{0}; i < nPoints; i++) {
    double E{0.1 + endpointEnergy * double(i) / double(nPoints - 1)};
    energyVec.push_back(E);
    decayRateVec.push_back(TritiumDecayRate(E, 0, 0, 0, endpointEnergy));
  }
  return pld(energyVec.begin(), energyVec.end(), decayRateVec.begin());
}

// Main function
int main(int argc, char* argv[]) {
  int opt{};
  std::string outputFile{" "};
  unsigned int nSims{100};
  while ((opt = getopt(argc, argv, "o:n:")) != -1) {
    switch (opt) {
      case 'o':
        outputFile = optarg;
        break;
      case 'n':
        nSims = std::stoi(optarg);
        break;
      case ':':
        std::cerr << "Option -" << optopt << " requires an argument."
                  << std::endl;
        return 1;
      case '?':
        std::cerr << "Unknown option: " << optopt << std::endl;
        return 1;
      default:
        std::cerr << "Usage: " << argv[0] << " [-o outputfile] [-n nSims]"
                  << std::endl;
        return 1;
    }
  }

  std::cout << "Output file: " << outputFile << std::endl;
  TString outputFileTStr{outputFile};
  TFile outFile(outputFileTStr, "RECREATE");
  TTree outTree("secElecTree", "Secondary electron production");

  double startKE{};
  double exitTime{};
  int nSecondaries{};
  const unsigned int nMaxElectrons{100};
  double secondaryKE[nMaxElectrons];
  outTree.Branch("startKE", &startKE, "startKE/D");
  outTree.Branch("exitTime", &exitTime, "exitTime/D");
  outTree.Branch("nSecondaries", &nSecondaries, "nSecondaries/I");
  outTree.Branch("secondaryKE", secondaryKE, "secondaryKE[nSecondaries]/D");

  const double molNumberDens{1e21};                       // m^-3
  const double THalfLife{12.32 * 365.25 * 24 * 60 * 60};  // seconds
  const double decayConst{log(2) / THalfLife};            // s^-1

  // Density decay rate
  const double densDecayRate{molNumberDens * 2 * decayConst};
  std::cout << "Density decay rate: " << densDecayRate * 1e-6
            << " cm^-3 s^-1 vs. literature value of 3.6e6 cm^-3 s^-1"
            << std::endl;

  // Define the WGTS
  const double WGTSLength{10};     // metres
  const double WGTSRadius{45e-3};  // metres
  const double WGTSField{3.6};     // Tesla
  auto field{new UniformField(WGTSField)};

  // Calc the step size we want here
  const double EMax{18.6e3};  // eV
  const double EMin{100.0};   // eV
  const double maxCycFreq{CalcCyclotronFreq(EMax, WGTSField)};
  std::cout << "Rough cyclotron frequency: " << maxCycFreq << " Hz"
            << std::endl;
  const double deltaT{1 / (10 * maxCycFreq)};  // seconds
  std::cout << "Time step: " << deltaT << " seconds" << std::endl;
  const double tau{2 * R_E / (3 * TMath::C())};

  // Set up random number stuff
  std::random_device rd{};
  std::mt19937 gen(rd());

  pld spectrum{GenerateSpectrum(10000)};

  for (unsigned int iEv{0}; iEv < nSims; iEv++) {
    double EGen{spectrum(gen)};
    startKE = EGen;
    std::cout << "Simulating event " << iEv + 1 << "\tE = " << EGen << " eV"
              << std::endl;

    nSecondaries = 0;
    for (unsigned int i{0}; i < nMaxElectrons; i++) {
      secondaryKE[i] = -1;
    }

    // Generate a random position in the cylinder
    std::uniform_real_distribution<double> uni1(0, 1);
    double zGen{0};
    double thetaPosGen{uni1(gen) * 2 * M_PI};
    double rGen{WGTSRadius * sqrt(uni1(gen))};

    TVector3 genPos(rGen * cos(thetaPosGen), rGen * sin(thetaPosGen), zGen);

    // Generate a random energy according to the tritium spectrum
    const double vGen{GetSpeedFromKE(EGen, ME)};
    // Generate isotropic velocity
    const double phiVelGen{uni1(gen) * 2 * M_PI};
    const double thetaVelGen{acos(2 * uni1(gen) - 1)};
    TVector3 genVel(vGen * sin(thetaVelGen) * cos(phiVelGen),
                    vGen * sin(thetaVelGen) * sin(phiVelGen),
                    vGen * cos(thetaVelGen));

    double tSim{0};              // seconds
    const double maxTSim{1e-5};  // seconds

    // Set up the Boris solver
    BorisSolver solver(field, -TMath::Qe(), ME, tau);

    std::vector<double> secondaryElecKE{};

    while (tSim < maxTSim && abs(genPos.Z()) < WGTSLength / 2 && EGen > 100) {
      // Set up the cross sections
      ElasticScatter elasticScatter(EGen);
      InelasticScatter inelasticScatter(EGen, CalcType::Kim1994, Species::H2);
      const double elasticXSec{elasticScatter.GetTotalXSec()};
      const double inelasticXSec{inelasticScatter.GetTotalXSec()};
      const double totalXSec{elasticXSec + inelasticXSec};
      const double meanFreeTime{1 / (totalXSec * molNumberDens * genVel.Mag())};
      // Draw the scattering time from an exponential distribution
      std::exponential_distribution<double> expDist(1 / meanFreeTime);
      const double tScatter{expDist(gen)};

      // Check if the scattering time is less than the time step
      if (tScatter < deltaT) {
        tSim += tScatter;
        auto [pos, vel] = solver.advance_step(tScatter, genPos, genVel);
        genPos = pos;
        genVel = vel;
        EGen = (sqrt(1 / (1 - pow(genVel.Mag() / TMath::C(), 2))) - 1) * ME_EV;

        // Check if this is an inelastic or elastic scattering event
        if (uni1(gen) < elasticXSec / totalXSec) {
          // Elastic scattering
          const double theta{elasticScatter.GetRandomScatteringAngle()};

          // Update the velocity
          genVel = elasticScatter.GetScatteredVector(genVel, EGen, theta);
        } else {
          // Inelastic scattering
          // Get the energy of the outgoing electron
          const double W{inelasticScatter.GetRandomW()};
          const double theta{inelasticScatter.GetRandomTheta(W)};

          secondaryElecKE.push_back(EGen - W - RYDBERG_EV);
          secondaryKE[nSecondaries] = EGen - W - RYDBERG_EV;
          nSecondaries++;
          if (nSecondaries >= nMaxElectrons) {
            std::cerr << "Too many secondary electrons produced" << std::endl;
            break;
          }

          // Update the energy and velocity
          genVel = inelasticScatter.GetScatteredVector(genVel, W, theta);
          EGen = W;
        }
      } else {
        auto [pos, vel] = solver.advance_step(deltaT, genPos, genVel);
        tSim += deltaT;
        genPos = pos;
        genVel = vel;
        EGen = (sqrt(1 / (1 - pow(genVel.Mag() / TMath::C(), 2))) - 1) * ME_EV;
      }
    }

    exitTime = tSim;
    std::cout << "Electron exited after " << tSim * 1e9 << " ns" << std::endl;
    std::cout << "Produced " << secondaryElecKE.size()
              << " secondary electrons with energies:" << std::endl;
    for (const auto& elecKE : secondaryElecKE) {
      std::cout << elecKE << " eV" << std::endl;
    }
    std::cout << "\n";

    outTree.Fill();
  }

  outFile.cd();
  outTree.Write();
  outFile.Close();

  delete field;
  return 0;
}