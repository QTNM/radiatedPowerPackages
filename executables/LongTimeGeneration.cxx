/*
  LongTimeGeneration.cxx

  Try and generate a nice track with >1ms total time
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
#include "TMath.h"
#include "TString.h"
#include "TSystem.h"
#include "TTree.h"
#include "TVector3.h"

using namespace rad;

using std::cout;
using std::endl;

int main(int argc, char *argv[]) {
  int opt{};
  std::string outputStemStr{" "};
  unsigned int nElectrons{2000};

  while ((opt = getopt(argc, argv, ":o:n:")) != -1) {
    switch (opt) {
      case 'o':
        outputStemStr = optarg;
        std::cout << "Output directory is " << outputStemStr << std::endl;
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

  std::cout << "Attempting to generate " << nElectrons << " electrons\n";

  // Define a magnetic trap. Start with a harmonic one
  const double coilRadius{0.03};                               // m
  const double trapDepth{4e-3};                                // T
  const double coilCurrent{2 * trapDepth * coilRadius / MU0};  // T
  const double bkg{0.7};                                       // Tesla
  auto field = new HarmonicField(coilRadius, coilCurrent, bkg);
  const double tritiumDensity{1e18};  // atoms m^-3
  const double zLimit{0.1};

  // Set up random number stuff
  std::random_device rd;
  std::mt19937 gen(rd());

  TString outputStem{outputStemStr};

  double simStepTime{5e-12};  // seconds
  const double tau{2 * R_E / (3 * TMath::C())};

  for (unsigned int i{0}; i < nElectrons; i++) {
    cout << "Electron " << i << endl;
    // Generate a random position in a cylinder
    const double rMax{0.015};  // m
    const double zMax{0.025};  // m

    double eKE{18.6e3};                      // eV
    double eSpeed{GetSpeedFromKE(eKE, ME)};  // m/s

    // Uniform distribution
    std::uniform_real_distribution<double> uni1(0, 1);
    double zGen{-zMax + uni1(gen) * 2 * zMax};
    double thetaPosGen{uni1(gen) * TMath::TwoPi()};
    double rGen{rMax * sqrt(uni1(gen))};

    TVector3 pos(rGen * cos(thetaPosGen), rGen * sin(thetaPosGen), zGen);

    // Generate a uniform velocity distribution
    double phiVelGen{uni1(gen) * TMath::TwoPi()};
    double thetaVelGen{uni1(gen) * TMath::Pi()};
    TVector3 vel(eSpeed * cos(phiVelGen) * sin(thetaVelGen),
                 eSpeed * sin(phiVelGen) * sin(thetaVelGen),
                 eSpeed * cos(thetaVelGen));

    // Generate a file and tree
    TString outputFile{outputStem + Form("/track%d.root", i)};
    TFile fout(outputFile, "recreate");
    auto tree = new TTree("tree", "tree");
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

    time = 0;
    xPos = pos.X();
    yPos = pos.Y();
    zPos = pos.Z();
    xVel = vel.X();
    yVel = vel.Y();
    zVel = vel.Z();

    BorisSolver solver(field, -TMath::Qe(), ME, tau);
    TVector3 acc{solver.acc(pos, vel)};
    xAcc = acc.X();
    yAcc = acc.Y();
    zAcc = acc.Z();
    tree->Fill();

    // Calculate next scattering time
    ElasticScatter scatEl(eKE);
    double elXSec{scatEl.GetTotalXSec()};
    InelasticScatter scatInel(eKE);
    double inelXSec{scatInel.GetTotalXSec()};
    double totalXSec{elXSec + inelXSec};
    // Calculate the mean free path
    const double lambdaStep{1 / (tritiumDensity * totalXSec)};
    std::exponential_distribution<double> pathDistStep(1 / lambdaStep);
    const double pathLenStep{pathDistStep(gen)};
    const double pathTimeStep{pathLenStep / vel.Mag()};
    double nextScatterTime{pathTimeStep};
    std::cout << "Scattering after " << pathTimeStep * 1e6 << " us\n";

    // Now propagate the electron
    unsigned long nSteps{1};
    while (abs(pos.Z()) < zLimit / 2) {
      // Step forward
      time = double(nSteps) * simStepTime;
      if (std::fmod(time, 5e-6) < simStepTime) {
        cout << "Simulated " << time * 1e6 << " us\n";
      }

      auto outputStep = solver.advance_step(simStepTime, pos, vel);
      pos = std::get<0>(outputStep);

      if (time < nextScatterTime) {
        vel = std::get<1>(outputStep);
        acc = solver.acc(pos, vel);
        // Update position and velocity vectors and write to tree
        xPos = pos.X();
        yPos = pos.Y();
        zPos = pos.Z();
        xVel = vel.X();
        yVel = vel.Y();
        zVel = vel.Z();
        xAcc = acc.X();
        yAcc = acc.Y();
        zAcc = acc.Z();
        tree->Fill();
      } else {
        // Time to do a scattering calculation
        double gamma{1 / sqrt(1 - pow(vel.Mag() / TMath::C(), 2))};
        eKE = (gamma - 1) * ME_EV;

        // Recalculate the cross sections based on the cross-sections
        ElasticScatter scatEl2(eKE);
        InelasticScatter scatInel2(eKE);
        elXSec = scatEl2.GetTotalXSec();
        inelXSec = scatInel2.GetTotalXSec();
        totalXSec = elXSec + inelXSec;

        // Figure out if this an elastic or inelastic scatter
        double scatAngle{0};
        if (uni1(gen) < elXSec / totalXSec) {
          // We have an elastic scatter
          cout << "Elastic scatter\n";
          // No energy loss so just get the scattering angle
          scatAngle = scatEl2.GetRandomScatteringAngle();
          vel = scatEl2.GetScatteredVector(vel, eKE, scatAngle);
        } else {
          // We have an inelastic scatter
          cout << "Inelastic scatter\n";
          // Get the energy of the secondary
          double wSample{scatInel2.GetRandomW()};
          double theta2Sample{scatInel2.GetRandomTheta(wSample)};

          // Now calculate the energy and the scattering angle of the primary
          scatAngle = scatInel2.GetPrimaryScatteredAngle(wSample, theta2Sample);
          eKE = scatInel2.GetPrimaryScatteredE(wSample, theta2Sample);
          vel = scatInel2.GetScatteredVector(vel, eKE, scatAngle);
        }
        cout << "Scattering angle = " << scatAngle * 180 / TMath::Pi()
             << " degrees\t New KE = " << eKE / 1e3 << " keV\n";

        // Now calculate the new acceleration
        acc = solver.acc(pos, vel);

        xPos = pos.X();
        yPos = pos.Y();
        zPos = pos.Z();
        xVel = vel.X();
        yVel = vel.Y();
        zVel = vel.Z();
        xAcc = acc.X();
        yAcc = acc.Y();
        zAcc = acc.Z();
        tree->Fill();

        // Get the next scatter time
        ElasticScatter scatElNext(eKE);
        double elXSecNext{scatElNext.GetTotalXSec()};
        InelasticScatter scatInelNext(eKE);
        double inelXSecNext{scatInelNext.GetTotalXSec()};
        double totalXSecNext{elXSecNext + inelXSecNext};
        double lambdaStepNext{1 / (tritiumDensity * totalXSecNext)};
        std::exponential_distribution<double> pathDistStepNext(1 /
                                                               lambdaStepNext);
        const double pathLenStepNext{pathDistStepNext(gen)};
        const double pathTimeStepNext{pathLenStepNext / vel.Mag()};
        nextScatterTime += pathTimeStepNext;
        std::cout << "Next scatter at " << nextScatterTime * 1e6 << " us\n";
      }

      nSteps++;
    }

    // We have now completed the motion
    cout << "Escaped after " << time * 1e6 << " us\n";

    fout.cd();
    tree->Write("", TObject::kOverwrite);
    delete tree;
    fout.Close();

    // We actually only want those electrons which have travelled for more than
    // 1ms
    if (time < 1e-3) {
      cout << "Too short, deleting file\n\n";
      gSystem->Exec("rm -f " + outputFile);
    }
  }

  return 0;
}