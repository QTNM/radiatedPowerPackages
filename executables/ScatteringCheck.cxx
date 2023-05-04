/*
  ScatteringCheck.cxx
*/

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

int main() {
  const TString outputFile{
      "/home/sjones/work/qtnm/outputs/ScatteringCheck/plots.root"};
  auto fout = std::make_unique<TFile>(outputFile, "recreate");

  std::random_device rd;
  std::mt19937 gen(rd());

  // First generate the spectrum distribution
  const int nSpecPnts{2000};
  pld spec{GenSpectrum(nSpecPnts)};

  // Number of electrons to generate
  const unsigned int nElectrons{100};
  TH1D hEngInit("hEngInit", "Initial energies; E [keV]; N_{e}", 50, 0, 18.6);
  SetHistAttr(hEngInit);
  TH1D hPathLen("hPathLen", "Distance between scatters; L [m]; N_{e}", 50, 0,
                4000);
  SetHistAttr(hPathLen);
  TH1D hPathTime("hPathTime", "Time between scatters; t [s]; N_{e}", 50, 0,
                 70e-6);
  SetHistAttr(hPathTime);

  TH1D hScatters("hScatters",
                 "Number of scatters (trapped electrons); N; N_{e}", 10, 0, 10);
  SetHistAttr(hScatters);
  TH1D hLifetime("hLifetime", "Lifetime (trapped electrons); N_{e}", 50, 0,
                 500e-6);
  SetHistAttr(hLifetime);
  TH1D hElastic("hElastic", "Is scatter elastic?; Elastic?; N_{scatters}", 2,
                -0.5, 1.5);
  SetHistAttr(hElastic);
  TH1D hELoss("hELoss",
              "Energy loss (inelastic scatters); E [eV]; N_{scatters}", 50, 0,
              100);
  SetHistAttr(hELoss);
  TH1D hScatEl(
      "hScatEl",
      "Scattering angle (elastic scatters); #theta [degrees]; N_{scatters}", 50,
      0, 5);
  SetHistAttr(hScatEl);
  TH1D hScatInel(
      "hScatInel",
      "Scattering angle (inelastic scatters); #theta [degrees]; N_{scatters}",
      50, 0, 5);
  SetHistAttr(hScatInel);
  // Define a magnetic trap. Start with a harmonic one
  const double coilRadius{0.03};                               // m
  const double trapDepth{4e-3};                                // T
  const double coilCurrent{2 * trapDepth * coilRadius / MU0};  // T
  const double bkg{0.72};                                      // Tesla
  auto field = new HarmonicField(coilRadius, coilCurrent, bkg);
  const double tritiumDensity{2e18};                           // atoms m^-3
  const double zMax{0.1};                                      // metres

  for (unsigned int i{0}; i < nElectrons; i++) {
    int nScatters{0};
    TVector3 pos(0, 0, 0);  // metres
    TVector3 vel{GenDecayVelocity(spec, gen)};
    double speed{vel.Mag()};
    // Recover the electron kinetic energy
    double gamma{1 / sqrt(1 - pow(vel.Mag() / TMath::C(), 2))};
    double ke{(gamma - 1) * ME_EV};
    cout << "Electron " << i + 1 << ": E_i = " << ke / 1e3 << " keV\n";
    hEngInit.Fill(ke / 1e3);

    bool isTrapped{true};
    const double keCutoff{1};       // eV
    const double timeCutoff{1e-3};  // s
    double totalTime{0};            // s
    const double tau{2 * R_E / (3 * TMath::C())};

    double scatAngle{0};
    while (isTrapped && ke > 1 && totalTime < timeCutoff) {
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
      hPathLen.Fill(pathLenStep);
      hPathTime.Fill(pathTimeStep);

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
        nScatters++;
        // Now figure out kinematics of the scatter
        // Calculate the current kinetic energy
        gamma = 1 / sqrt(1 - pow(vel.Mag() / TMath::C(), 2));
        ke = (gamma - 1) * ME_EV;
        speed = vel.Mag();

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
        if (uni(gen) < elXSec / totalXSec) {
          // We have an elastic scatter
          cout << "Elastic scatter\n";
          // No energy loss so just get the scattering angle
          scatAngle = scatEl2.GetRandomScatteringAngle();
          hElastic.Fill(1);
          hScatEl.Fill(scatAngle);
        } else {
          // We have an inelastic scatter
          cout << "Inelastic scatter\n";
          // Get the energy of the secondary
          double wSample{scatInel2.GetRandomW()};
          double theta2Sample{scatInel2.GetRandomTheta(wSample)};
          // Now calculate the energy and the scattering angle of the primary
          scatAngle = scatInel2.GetPrimaryScatteredAngle(wSample, theta2Sample);
          hELoss.Fill(ke -
                      scatInel2.GetPrimaryScatteredE(wSample, theta2Sample));
          ke = scatInel2.GetPrimaryScatteredE(wSample, theta2Sample);
          hElastic.Fill(0);
          hScatInel.Fill(scatAngle);
        }
        speed = GetSpeedFromKE(ke, ME);

        cout << "Scattering angle = " << scatAngle * 180 / TMath::Pi()
             << " degrees\tNew KE = " << ke / 1e3 << " keV\n";

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
      }
    }

    // We have now completed the motion
    cout << "Escaped after " << totalTime * 1e6 << " us and " << nScatters
         << " scatters\n\n";

    if (nScatters >= 1) {
      hScatters.Fill(nScatters);
      hLifetime.Fill(totalTime);
    }
  }
  delete field;

  fout->cd();
  hEngInit.Write();
  hPathLen.Write();
  hPathTime.Write();

  hScatters.Write();
  hLifetime.Write();
  hElastic.Write();
  hScatEl.Write();
  hScatInel.Write();
  hELoss.Write();

  fout->Close();
  return 0;
}