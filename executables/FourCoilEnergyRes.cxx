/*
  FourCoilEnergyRes.cxx

  Generate tracks in multiple different configurations
*/

#include <getopt.h>
#include <unistd.h>

#include <boost/lexical_cast.hpp>
#include <boost/uuid/uuid.hpp>
#include <boost/uuid/uuid_generators.hpp>
#include <boost/uuid/uuid_io.hpp>
#include <ctime>
#include <iostream>
#include <random>
#include <string>
#include <tuple>

#include "BasicFunctions/BasicFunctions.h"
#include "BasicFunctions/Constants.h"
#include "ElectronDynamics/BorisSolver.h"
#include "ElectronDynamics/QTNMFields.h"
#include "SignalProcessing/LocalOscillator.h"
#include "SignalProcessing/Signal.h"
#include "TF1.h"
#include "TFile.h"
#include "TGraph.h"
#include "TString.h"
#include "TSystem.h"
#include "TTree.h"
#include "Waveguides/CircularWaveguide.h"

using namespace rad;

std::string make_uuid() {
  return boost::lexical_cast<std::string>((boost::uuids::random_generator())());
}

double UniformBkg(double *x, double *par) {
  // Trap coil definitions
  const double rCoil{15e-3};                                            // m
  const double outerCoilOffset{25e-3};                                  // m
  const double outerCoilTrapDepth{5.5e-3};                              // T
  const double outerCoilCurrent{2 * outerCoilTrapDepth * rCoil / MU0};  // A
  const double innerCoilOffset{19e-3};                                  // m
  const double innerCoilTrapDepth{3.1e-3};                              // T
  const double innerCoilCurrent{2 * innerCoilTrapDepth * rCoil / MU0};  // A

  // Set up uniform background field first
  const double bkgField{0.7};  // Tesla
  FourCoilBathtub uniformBkg(rCoil, outerCoilCurrent, outerCoilOffset, rCoil,
                             -innerCoilCurrent, innerCoilOffset, bkgField);

  double zPos{x[0] / 1e3};
  TVector3 pos(0, 0, zPos);
  return uniformBkg.evaluate_field_magnitude(pos);
}

double HelmholtzBkg(double *x, double *par) {
  // Trap coil definitions
  const double rCoil{15e-3};                                            // m
  const double outerCoilOffset{25e-3};                                  // m
  const double outerCoilTrapDepth{5.5e-3};                              // T
  const double outerCoilCurrent{2 * outerCoilTrapDepth * rCoil / MU0};  // A
  const double innerCoilOffset{19e-3};                                  // m
  const double innerCoilTrapDepth{3.1e-3};                              // T
  const double innerCoilCurrent{2 * innerCoilTrapDepth * rCoil / MU0};  // A

  // Helmholtz definitions
  const double n{430};
  const double lSolenoid{5e-2};
  const double rSolenoid{200e-3};
  const double bField{0.4912};
  const double iSolenoid{2 * bField *
                         sqrt(pow(lSolenoid / 2, 2) + rSolenoid * rSolenoid) /
                         (MU0 * n * lSolenoid)};  // Amps

  HelmholtzField bkg(rSolenoid, lSolenoid, iSolenoid, n);
  FourCoilField trap(rCoil, outerCoilCurrent, outerCoilOffset, rCoil,
                     -innerCoilCurrent, innerCoilOffset);
  FourCoilHelmholtz combinedField(trap, bkg);

  double zPos{x[0] / 1e3};
  TVector3 pos(0, 0, zPos);
  return combinedField.evaluate_field_magnitude(pos);
}

double BathtubAnalytical(double *x, double *par) {
  double z{x[0] / 1e3};
  double b0{par[0]};
  double l0{par[1]};
  double l1{par[2]};
  if (z < -l1 / 2) {
    return b0 * (1 + pow(z + l1 / 2, 2) / (l0 * l0));
  } else if (z >= -l1 / 2 && z <= l1 / 2) {
    return b0;
  } else {
    return b0 * (1 + pow(z - l1 / 2, 2) / (l0 * l0));
  }
}

double CalcAvgBathtubFreq(double ke, double b0, double l0, double l1,
                          double theta) {
  const double gamma{(ke + ME_EV) / ME_EV};
  const double prefactor{TMath::Qe() * b0 / (gamma * ME * 2 * M_PI)};
  const double zMax{l0 / tan(theta)};
  return prefactor * (1 + (zMax * zMax / (2 * l0 * l0)) *
                              pow(1 + l1 * tan(theta) / (M_PI * l0), -1));
}

double CalcBathtubAxFreq(double v0, double b0, double l0, double l1,
                         double theta) {
  double fa{v0 * sin(theta) / (2 * M_PI * l0)};
  return fa * pow(1 + l1 * tan(theta) / (l0 * M_PI), -1);
}

int main(int argc, char *argv[]) {
  // Start a clock
  const clock_t startGenClock{clock()};

  int opt{};
  std::string outputStemStr{" "};
  unsigned int nEvents{1000};
  bool uniformField{false};  // Boolean controlling the field used
  double maxSimTime{10e-6};  // Maximum time in seconds to simulate electron for

  while ((opt = getopt(argc, argv, ":o:n:t:f")) != -1) {
    switch (opt) {
      case 'o':
        outputStemStr = optarg;
        break;

      case 'n':
        nEvents = std::stoi(optarg);
        break;

      case 't':
        maxSimTime = std::stod(optarg);
        break;

      case 'f':
        uniformField = true;
        break;

      case ':':
        std::cout << "Option needs a value\n";
        break;

      case '?':
        std::cout << "Unknown option: " << optopt << std::endl;
        break;
    }
  }

  // Do some output checking and print out the basic simulation parameters
  std::cout << "Attempting to generate " << nEvents << " events/electrons.\n";
  std::cout << "Outputting files in directory " << outputStemStr << std::endl;
  std::cout << "Maximum simulation time is " << maxSimTime * 1e6 << " us\n";
  if (uniformField) {
    std::cout << "Using uniform background field\n";
  } else {
    std::cout << "Using Helmholtz background field\n";
  }

  // Create output file for plots
  std::string randomPlotFileName{make_uuid()};
  TString outputStem{outputStemStr};
  TString outputFile{outputStem +
                     Form("/plots_%s.root", randomPlotFileName.data())};
  TFile fout(outputFile, "recreate");

  TF1 fUniformBkg("fUniformBkg", UniformBkg, -50, 50, 0);
  fUniformBkg.SetNpx(1000);
  fUniformBkg.SetLineWidth(3);
  fUniformBkg.SetTitle("Uniform background; z [mm]; |B| [T]");

  TF1 fHelmholtzBkg("fHelmholtzBkg", HelmholtzBkg, -50, 50, 0);
  fHelmholtzBkg.SetNpx(1000);
  fHelmholtzBkg.SetLineWidth(3);
  fHelmholtzBkg.SetTitle("Helmholtz background; z [mm]; |B| [T]");
  fHelmholtzBkg.SetLineColor(kBlue);

  const double b0{0.70002};  // T
  const double l0{0.21};     // m
  const double l1{0.02};     // m
  TF1 fAna("fAna", BathtubAnalytical, -30, 30, 3);
  fAna.SetParameter(0, b0);
  fAna.SetParameter(1, l0);
  fAna.SetParameter(2, l1);

  fout.cd();
  fUniformBkg.Write();
  fHelmholtzBkg.Write();
  fAna.Write();

  // Define the trap field quantities
  const double rCoil{15e-3};                                            // m
  const double outerCoilOffset{25e-3};                                  // m
  const double outerCoilTrapDepth{5.5e-3};                              // T
  const double outerCoilCurrent{2 * outerCoilTrapDepth * rCoil / MU0};  // A
  const double innerCoilOffset{19e-3};                                  // m
  const double innerCoilTrapDepth{3.1e-3};                              // T
  const double innerCoilCurrent{2 * innerCoilTrapDepth * rCoil / MU0};  // A

  // Uniform field generation
  const double bkgField{0.7};  // Tesla

  // Generate the field depending on the boolean
  BaseField *field = 0;
  if (uniformField) {
    field = new FourCoilBathtub(rCoil, outerCoilCurrent, outerCoilOffset, rCoil,
                                -innerCoilCurrent, innerCoilOffset, 0.7);
  } else {
    // Helmholtz values
    const double n{430};
    const double lSolenoid{5e-2};
    const double rSolenoid{200e-3};
    const double bField{0.4912};
    const double iSolenoid{2 * bField *
                           sqrt(pow(lSolenoid / 2, 2) + rSolenoid * rSolenoid) /
                           (MU0 * n * lSolenoid)};  // Amps
    FourCoilField trapField(rCoil, outerCoilCurrent, outerCoilOffset, rCoil,
                            -innerCoilCurrent, innerCoilOffset);
    HelmholtzField bkg(rSolenoid, lSolenoid, iSolenoid, n);
    field = new FourCoilHelmholtz(trapField, bkg);
  }

  // Set up a random number generator
  std::random_device rd;
  std::mt19937 gen(rd());

  // About 15 steps per orbit
  const double simStepSize{3.5e-12};  // seconds
  // Sim time should give us 50 kHz frequency bin width

  // Set up constant electron kinematics
  const double eKE{18.6e3};  // eV
  const double eSpeed{GetSpeedFromKE(eKE, ME)};
  const double tau{2 * R_E / (3 * TMath::C())};

  // Generation specifics
  const double rGenMax{6e-3};  // metres
  const double zGen{0};

  auto trResult = new TTree("trResult", "trResult");
  const double minPeakSize{1e-13};  // Min power to record
  const int nMaxPeaks{100};
  unsigned int nPeaks{};
  double peakPos[nMaxPeaks];
  double peakHeight[nMaxPeaks];
  double theta{};
  double xStart{}, yStart{};
  double loFreq{};
  double calcAxFreq{};
  trResult->Branch("theta", &theta);
  trResult->Branch("xStart", &xStart);
  trResult->Branch("yStart", &yStart);
  trResult->Branch("loFreq", &loFreq);
  trResult->Branch("calcAxFreq", &calcAxFreq);
  trResult->Branch("nPeaks", &nPeaks, "nPeaks/i");
  trResult->Branch("peakPos", peakPos, "peakPos[nPeaks]/D");
  trResult->Branch("peakHeight", peakHeight, "peakHeight[nPeaks]/D");

  // Loop over the number of events, generating one electron at a time
  for (uint iEv{0}; iEv < nEvents; iEv++) {
    const clock_t eventStartClock{clock()};

    // Create a uniform distribution
    std::uniform_real_distribution<double> uni1(0, 1);
    // Generate the initial electron position
    const double rGen{rGenMax * sqrt(uni1(gen))};
    const double thetaPosGen{uni1(gen) * 2 * M_PI};
    TVector3 x0(rGen * cos(thetaPosGen), rGen * sin(thetaPosGen), zGen);

    // Generate the initial electron velocity
    const double phiVelGen{uni1(gen) * 2 * M_PI};
    const double thetaVelGen{uni1(gen) * M_PI};
    TVector3 v0(eSpeed * cos(phiVelGen) * sin(thetaVelGen),
                eSpeed * sin(phiVelGen) * sin(thetaVelGen),
                eSpeed * cos(thetaVelGen));
    const double thetaBotRad{abs(atan(v0.Perp() / v0.Z()))};
    const double thetaBotDeg{thetaBotRad * 180 / M_PI};
    std::cout << "Electron " << iEv << ": theta = " << thetaBotDeg
              << " degrees\n";

    // Generate a file and accompanying tree
    std::string randomTrackFileName{make_uuid()};
    TString trackFile{outputStem +
                      Form("/track_%s.root", randomTrackFileName.data())};
    TFile fTrack(trackFile, "recreate");
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
    xPos = x0.X();
    yPos = x0.Y();
    zPos = x0.Z();
    xVel = v0.X();
    yVel = v0.Y();
    zVel = v0.Z();
    BorisSolver solver(field, -TMath::Qe(), ME, tau);
    TVector3 acc{solver.acc(x0, v0)};
    xAcc = acc.X();
    yAcc = acc.Y();
    zAcc = acc.Z();
    tree->Fill();

    // Now propagate the electron
    TVector3 vel{v0};
    TVector3 pos{x0};
    unsigned long nSteps{1};
    bool isTrapped{true};
    while (time <= maxSimTime && isTrapped) {
      time = double(nSteps) * simStepSize;
      // std::cout << nSteps << " steps\t" << time << std::endl;
      if (std::fmod(time, 1e-6) < simStepSize) {
        std::cout << "Simulated " << time * 1e6 << " us\n";
      }

      // Advance the electron
      auto outputStep = solver.advance_step(simStepSize, pos, vel);
      pos = std::get<0>(outputStep);
      if (abs(pos.Z()) > outerCoilOffset) {
        std::cout << "Electron escaped after " << time * 1e6 << " us\n";
        isTrapped = false;
        break;
      }

      vel = std::get<1>(outputStep);
      acc = solver.acc(pos, vel);

      // Write to tree
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
      nSteps++;
    }

    // Write the tree to file
    fTrack.cd();
    tree->Write("", TObject::kOverwrite);
    delete tree;
    fTrack.Close();

    // If the electron is trapped, run the signal processing
    if (isTrapped) {
      std::cout << "Running signal processing\n";
      // Generate the waveguide we want
      const double wgRadius{7.14e-3};         // metres
      const double wgLength{0.08};            // metres
      TVector3 probePos(0, 0, wgLength / 2);  // At the end of the guide
      auto wg = new CircularWaveguide(wgRadius, wgLength, probePos);

      const double sampleRate{1.5e9};  // Hertz
      // Set up the local oscillator
      const double avgCycFreq{CalcAvgBathtubFreq(eKE, b0, l0, l1, thetaBotRad)};
      calcAxFreq = CalcBathtubAxFreq(eSpeed, b0, l0, l1, thetaBotRad);

      std::cout << "Average cyclotron frequency = " << avgCycFreq / 1e9
                << " GHz\tAxial frequency = " << calcAxFreq / 1e6 << " MHz\n";
      loFreq = avgCycFreq - sampleRate / 4;
      LocalOscillator lo(2 * M_PI * loFreq);
      // Actually generate the signal
      Signal sig(trackFile, wg, lo, sampleRate);
      auto grV{sig.GetVIPowerPeriodogram(1)};
      fout.cd();
      grV->Write(Form("grV_%d", iEv));

      // Set result tree
      theta = thetaBotRad;
      xStart = x0.X();
      yStart = x0.Y();

      // Reset the arrays
      for (uint iArr{0}; iArr < nMaxPeaks; iArr++) {
        peakHeight[iArr] = -1;
        peakPos[iArr] = -1;
      }

      // Extract peak information
      nPeaks = 0;
      for (int iN{1}; iN < grV->GetN() - 1; iN++) {
        if (grV->GetPointY(iN) < minPeakSize)
          continue;
        else if (grV->GetPointY(iN) >= minPeakSize &&
                 grV->GetPointY(iN) > grV->GetPointY(iN - 1) &&
                 grV->GetPointY(iN) > grV->GetPointY(iN + 1)) {
          double peakSum{grV->GetPointY(iN) + grV->GetPointY(iN + 1) +
                         grV->GetPointY(iN - 1)};
          peakHeight[nPeaks] = peakSum;
          peakPos[nPeaks] = grV->GetPointX(iN);
          nPeaks++;
        }
      }

      trResult->Fill();

      delete grV;
      delete wg;
    }

    // Delete the track file
    gSystem->Exec("rm -f " + trackFile);

    // Measure how long we've been generating for
    const clock_t eventEndClock{clock()};
    double eventTime{double(eventEndClock - eventStartClock) / CLOCKS_PER_SEC};
    std::cout << "Took " << eventTime << " seconds to generate event\n\n";
  }

  fout.cd();
  trResult->Write("", TObject::kOverwrite);

  fout.Close();
  delete field;

  const clock_t totalClock{clock()};
  const double totalTime{double(totalClock - startGenClock) / CLOCKS_PER_SEC};
  std::cout << "Total execution time is " << totalTime << " seconds\n";
  return 0;
}