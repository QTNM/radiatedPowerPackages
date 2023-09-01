/// crudeEnergyResolution.cxx
/// Program to measure the energy resolution of endpoint electrons for a given
/// trap geometry Electrons are generated isotropically throughout the trap and
/// the trapping condition checked

#include "Antennas/HalfWaveDipole.h"
#include "BasicFunctions/Constants.h"
#include "ElectronDynamics/BorisSolver.h"
#include "ElectronDynamics/QTNMFields.h"
#include "SignalProcessing/InducedVoltage.h"
#include "SignalProcessing/Signal.h"

// STL includes
#include <getopt.h>
#include <unistd.h>

#include <cmath>
#include <ctime>
#include <filesystem>
#include <iostream>
#include <string>
#include <tuple>

// ROOT includes
#include "TFile.h"
#include "TMath.h"
#include "TRandom3.h"
#include "TSystem.h"
#include "TTree.h"
#include "TVector3.h"

using namespace rad;

void PrintHelp() {
  std::cout
      << "--radius <r>:      Set the coil radius in metres\n"
         "--length <l>:      Set the trap length in metres\n"
         "--number <n>:      Number of electrons to simulate\n"
         "--inhomAx:     Fractional inhomogeneity in the axial direction\n"
         "--inhomRad:    Fractional inhomogeneity in the radial direction\n"
         "--outputDir <d>:   Output directory to write to\n"
         "--keepTracks <k>:  Sets flag to keep electron trajectories\n"
         "--help <h>:        Prints this help message\n";
  exit(1);
}

int main(int argc, char* argv[]) {
  int opt;

  double trapLength = 0.3;  // m
  double RCoil = 0.025;     // m
  int nElectrons = 1000;
  std::string outputDir = " ";
  bool keepTracks = false;
  double inhomAx = 0.0;
  double inhomRad = 0.0;

  const option long_opts[] = {{"outputDir", required_argument, nullptr, 'd'},
                              {"number", required_argument, nullptr, 'n'},
                              {"radius", required_argument, nullptr, 'r'},
                              {"length", required_argument, nullptr, 'l'},
                              {"inhomAx", required_argument, nullptr, 'z'},
                              {"inhomRad", required_argument, nullptr, 'x'},
                              {"keepTracks", no_argument, nullptr, 'k'},
                              {"help", no_argument, nullptr, 'h'},
                              {nullptr, no_argument, nullptr, 0}};

  while ((opt = getopt_long(argc, argv, ":d:n:r:l:kh", long_opts, nullptr)) !=
         -1) {
    switch (opt) {
      case 'h':
        PrintHelp();
      case 'd':
        outputDir = optarg;
        break;
      case 'n':
        nElectrons = atoi(optarg);
        break;
      case 'r':
        RCoil = atof(optarg);
        break;
      case 'l':
        trapLength = atof(optarg);
        break;
      case 'z':
        inhomAx = atof(optarg);
        break;
      case 'x':
        inhomRad = atof(optarg);
        break;
      case 'k':
        keepTracks = true;
        std::cout
            << "WARNING!: Choosing to keep the track output. The full tracks "
               "take up about 1GB each so make sure you have enough space."
            << std::endl;
        break;
      case ':':
        std::cout << "Option needs a value" << std::endl;
        break;
      case '?':
        std::cout << "Unknown option: " << optopt << std::endl;
        break;
    }
  }

  // Check parameter values
  if (nElectrons <= 0) {
    std::cout << "Invalid simulated number of electrons" << std::endl;
    exit(1);
  }
  if (outputDir == " ") {
    std::cout << "Must specify directory file with -d" << std::endl;
    exit(1);
  }
  if (RCoil <= 0) {
    std::cout << "Invalid coil radius provided" << std::endl;
    exit(1);
  }
  if (trapLength <= 0) {
    std::cout << "Invalid trap length provided" << std::endl;
    exit(1);
  }

  std::cout << "Output directory is " << outputDir << std::endl;
  std::cout << "Simulating " << nElectrons << " electrons" << std::endl;
  std::cout << "Coil radius is " << RCoil << " m" << std::endl;
  std::cout << "Chosen trap length is " << trapLength << " m" << std::endl;
  std::cout << "Axial inhomogeneity is " << inhomAx << std::endl;
  std::cout << "Radial inhomogeneity is " << inhomRad << std::endl;

  // RNG
  TRandom3* thisRand = new TRandom3(0);

  // Simulation parameters
  const double timeStepSize = 3.7e-12;  // seconds
  const double maxSimTime = 55e-6;      // seconds
  const double RGen = 0.02;             // m
  const double trapDepth = 0.0049;      // Tesla
  const double centralField = 1.0;      // Tesla

  InhomogeneousBackgroundField* bkg = new InhomogeneousBackgroundField(
      centralField, inhomAx, trapLength / 2, inhomRad, RCoil);
  const double trapFieldOffset =
      centralField -
      bkg->evaluate_field_at_point(TVector3(0, 0, trapLength / 2)).Mag();

  const double ICoil =
      2.0 * (trapDepth + trapFieldOffset) * RCoil / MU0;  // Amps
  delete bkg;

  // Generate the bathtub field
  InhomogeneousBathtubField* bathtubField = new InhomogeneousBathtubField(
      RCoil, ICoil, trapLength / 2, centralField, inhomAx, inhomRad);
  // Antenna specifications
  const double antennaRadius = 0.03;
  const double antennaAngle1 = 0 * TMath::Pi() / 180;
  const double antennaAngle2 = 90 * TMath::Pi() / 180;
  TVector3 antennaPoint1(antennaRadius * TMath::Cos(antennaAngle1),
                         antennaRadius * TMath::Sin(antennaAngle1), 0.0);
  TVector3 antennaDirZ1(-1 * TMath::Sin(antennaAngle1),
                        TMath::Cos(antennaAngle1), 0.0);
  TVector3 antennaDirX1(TMath::Cos(antennaAngle1), TMath::Sin(antennaAngle1),
                        0.0);
  HalfWaveDipole* antenna1 =
      new HalfWaveDipole(antennaPoint1, antennaDirX1, antennaDirZ1, 27.01e9);
  const double antennaLowerBandwidth = 26.95e9;
  const double antennaUpperBandwidth = 27.05e9;
  antenna1->SetBandwidth(antennaLowerBandwidth, antennaUpperBandwidth);

  // Electron dynamics
  const double TElec = 18600;  // eV
  const double gamma = TElec * TMath::Qe() / (ME * TMath::C() * TMath::C()) + 1;
  const double betaSq = 1 - 1 / pow(gamma, 2);
  const double initialSpeed = sqrt(betaSq) * TMath::C();
  const double tau = 2 * R_E / (3 * TMath::C());

  // Parameters for signal processing
  const double tAcq = maxSimTime - 1e-6;  // seconds
  const double loFreq = 26.75e9;          // Hz
  const double loadResistance = 70;       // Ohms
  const double sampleRate = 750e6;        // Hz
  LocalOscillator lo(2 * TMath::Pi() * loFreq);

  // Create the output directory if it doesn't exist
  bool directoryExists = std::filesystem::is_directory(outputDir);
  if (!directoryExists) {
    // Attempt to create the directory
    bool directoryCreated = std::filesystem::create_directories(outputDir);
    if (!directoryCreated) {
      std::cout << "Failed to create output directory" << std::endl;
      exit(2);
    }
  }

  TFile* fout = new TFile(Form("%s/output.root", outputDir.data()), "RECREATE");

  TTree* startTree = new TTree("startTree", "startTree");
  double initialAngle;
  double pitchAngle;
  double BMean;
  double xpos, ypos, zpos;
  double xvel, yvel, zvel;
  int isTrapped_val;
  startTree->Branch("initialAngle", &initialAngle, "initialAngle/D");
  startTree->Branch("pitchAngle", &pitchAngle, "pitchAngle/D");
  startTree->Branch("BMean", &BMean, "BMean/D");
  startTree->Branch("xpos", &xpos, "xpos/D");
  startTree->Branch("ypos", &ypos, "ypos/D");
  startTree->Branch("zpos", &zpos, "zpos/D");
  startTree->Branch("xvel", &xvel, "xvel/D");
  startTree->Branch("yvel", &yvel, "yvel/D");
  startTree->Branch("zvel", &zvel, "zvel/D");
  startTree->Branch("isTrapped", &isTrapped_val, "isTrapped/I");

  TH1D* hInitialAngle = new TH1D(
      "hInitialAngle", "Initial angles; Angle [degrees]; N", 180, 0, 90);
  SetHistAttr(hInitialAngle);
  TH1D* hPitchAngle =
      new TH1D("hPitchAngle", "Pitch angles; Angle [degrees]; N", 180, 0, 90);
  SetHistAttr(hPitchAngle);
  TH1D* hInitialAngleAcc = new TH1D(
      "hInitialAngleAcc", "Initial angles; Angle [degrees]; N", 90, 0, 90);
  SetHistAttr(hInitialAngleAcc);
  TH1D* hPitchAngleAcc =
      new TH1D("hPitchAngleAcc", "Pitch angles; Angle [degrees]; N", 90, 0, 90);
  SetHistAttr(hPitchAngleAcc);
  TH1D* hRPos = new TH1D("hRPos", "Radial position; R [m]; N", 50, 0, RGen);
  SetHistAttr(hRPos);
  TH1D* hRPosAcc =
      new TH1D("hRPosAcc", "Radial position; R [m]; N", 50, 0, RGen);
  SetHistAttr(hRPosAcc);
  TH1D* hZPos = new TH1D("hZPos", "Axial position; z [m]; N", 50,
                         -trapLength / 2, trapLength / 2);
  SetHistAttr(hZPos);
  TH1D* hZPosAcc = new TH1D("hZPosAcc", "Axial position; z [m]; N", 50,
                            -trapLength / 2, trapLength / 2);
  SetHistAttr(hZPosAcc);
  TH2D* h2RZPos = new TH2D("h2RZPos", "Position; z [m]; R [m]; N", 8,
                           -trapLength / 2, trapLength / 2, 8, 0, RGen);
  SetHistAttr(h2RZPos);
  TH2D* h2RZPosAcc = new TH2D("h2RZPosAcc", "Position; z [m]; R [m]; N", 8,
                              -trapLength / 2, trapLength / 2, 8, 0, RGen);
  SetHistAttr(h2RZPosAcc);

  // Generate the electrons
  for (int n = 0; n < nElectrons; n++) {
    std::cout << "Electron number " << n << std::endl;

    const clock_t begin_time = clock();

    // Generate random electron directions
    double phiVelGen = thisRand->Uniform() * 2 * TMath::Pi();
    double thetaVelGen = thisRand->Uniform() * TMath::Pi();
    // Generate random uniform volume distribution
    double radialPosGen = sqrt(thisRand->Uniform()) * RGen;
    double thetaPosGen = thisRand->Uniform() * 2 * TMath::Pi();
    double zPosGen = -trapLength / 2 + (thisRand->Uniform() * trapLength);

    TVector3 posVec(radialPosGen * cos(thetaPosGen),
                    radialPosGen * sin(thetaPosGen), zPosGen);
    TVector3 velVec(initialSpeed * cos(phiVelGen) * sin(thetaVelGen),
                    initialSpeed * sin(phiVelGen) * sin(thetaVelGen),
                    initialSpeed * cos(thetaVelGen));

    BMean = 0.0;
    int nRecordedSteps = 0;

    xpos = posVec.X();
    ypos = posVec.Y();
    zpos = posVec.Z();
    xvel = velVec.X();
    yvel = velVec.Y();
    zvel = velVec.Z();

    double RVel = sqrt(velVec.X() * velVec.X() + velVec.Y() * velVec.Y());
    std::cout << RVel << ", " << velVec.Z() << std::endl;
    initialAngle = abs(atan(RVel / velVec.Z()));
    pitchAngle = abs(atan(RVel / velVec.Z()));  // Initial value

    std::cout << "Initial angle (degrees), rPos, zPos = "
              << (initialAngle * 180 / TMath::Pi()) << ", " << radialPosGen
              << " m, " << zPosGen << " m" << std::endl;

    // Set up the solver
    BorisSolver solver(bathtubField, -TMath::Qe(), ME, tau);
    // Calculate the number of time steps
    int nTimeSteps = maxSimTime / timeStepSize;

    bool isTrapped = true;
    isTrapped_val = 1;

    hZPos->Fill(posVec.Z());
    hRPos->Fill(radialPosGen);
    h2RZPos->Fill(posVec.Z(), radialPosGen);
    hInitialAngle->Fill(initialAngle * 180 / TMath::Pi());

    TFile* fElec =
        new TFile(Form("%s/track%d.root", outputDir.data(), n), "RECREATE");

    // Create a new TTree representing the trajectory of this electron
    TTree* tree = new TTree("tree", "tree");
    double time;
    double xPos, yPos, zPos;
    double xVel, yVel, zVel;
    double xAcc, yAcc, zAcc;
    tree->Branch("time", &time, "time/D");
    tree->Branch("xPos", &xPos, "xPos/D");
    tree->Branch("yPos", &yPos, "yPos/D");
    tree->Branch("zPos", &zPos, "zPos/D");
    tree->Branch("xVel", &xVel, "xVel/D");
    tree->Branch("yVel", &yVel, "yVel/D");
    tree->Branch("zVel", &zVel, "zVel/D");
    tree->Branch("xAcc", &xAcc, "xAcc/D");
    tree->Branch("yAcc", &yAcc, "yAcc/D");
    tree->Branch("zAcc", &zAcc, "zAcc/D");

    // Set the initial state
    TVector3 eAcc = solver.acc(posVec, velVec);
    time = 0;
    xPos = posVec.X();
    yPos = posVec.Y();
    zPos = posVec.Z();
    xVel = velVec.X();
    yVel = velVec.Y();
    zVel = velVec.Z();
    xAcc = eAcc.X();
    yAcc = eAcc.Y();
    zAcc = eAcc.Z();

    tree->Fill();

    for (int iStep = 0; iStep < nTimeSteps; iStep++) {
      std::tuple<TVector3, TVector3> outputStep =
          solver.advance_step(timeStepSize, posVec, velVec);
      posVec = std::get<0>(outputStep);
      velVec = std::get<1>(outputStep);
      eAcc = solver.acc(posVec, velVec);

      time = double(iStep + 1) * timeStepSize;
      if (std::fmod(time, 5e-6) < timeStepSize)
        std::cout << time << " seconds generated" << std::endl;

      xPos = posVec.X();
      yPos = posVec.Y();
      zPos = posVec.Z();
      xVel = velVec.X();
      yVel = velVec.Y();
      zVel = velVec.Z();
      xAcc = eAcc.X();
      yAcc = eAcc.Y();
      zAcc = eAcc.Z();

      tree->Fill();

      double thisPitchAngle = abs(atan(RVel / velVec.Z()));
      if (thisPitchAngle < pitchAngle) pitchAngle = thisPitchAngle;

      // Check if the electron has escaped the trap
      if (abs(posVec.Z()) > trapLength / 2) {
        isTrapped = false;
        isTrapped_val = 0;
        std::cout << "Was not trapped" << std::endl;
        break;
      }

      double BField = bathtubField->evaluate_field_at_point(posVec).Mag();
      BMean += BField;
      nRecordedSteps++;
    }

    BMean /= double(nRecordedSteps);

    hPitchAngle->Fill(pitchAngle * 180 / TMath::Pi());

    fElec->cd();
    tree->Write();
    delete tree;
    fElec->Close();
    delete fElec;

    if (isTrapped) {
      std::cout << "Was trapped" << std::endl;
      hZPosAcc->Fill(posVec.Z());
      hRPosAcc->Fill(radialPosGen);
      h2RZPosAcc->Fill(posVec.Z(), radialPosGen);
      hInitialAngleAcc->Fill(initialAngle * 180 / TMath::Pi());
      hPitchAngleAcc->Fill(pitchAngle * 180 / TMath::Pi());

      // If the electron was trapped we can do some signal processing
      TString trackFile{Form("%s/track%d.root", outputDir.data(), n)};
      // Create the Signal with no noise
      Signal sig(trackFile, antenna1, lo, sampleRate, {});

      // Get the spectrogram and time domain signal
      TGraph* grVI = sig.GetVITimeDomain();
      grVI->SetTitle(Form(
          "Length = %.2f m: #theta = %.2f deg, r_{i} = %.2f m, z_{i} = %.2f m",
          trapLength, (pitchAngle * 180 / TMath::Pi()), radialPosGen, zPosGen));
      TGraph* grVIPgram = sig.GetVIPowerPeriodogram(loadResistance);
      grVIPgram->SetTitle(Form(
          "Length = %.2f m: #theta = %.2f deg, r_{i} = %.2f m, z_{i} = %.2f m",
          trapLength, (pitchAngle * 180 / TMath::Pi()), radialPosGen, zPosGen));

      fout->cd();
      grVI->Write(Form("grVI%d", n));
      grVIPgram->Write(Form("grVIPgram%d", n));
      delete grVI;
      delete grVIPgram;
    }  // Signal processing if electron is trapped

    const clock_t end_time = clock();
    std::cout << "Simulation time was "
              << float(end_time - begin_time) / CLOCKS_PER_SEC << " seconds"
              << std::endl;

    startTree->Fill();

    // If chosen, delete the file containing the track
    // This is necessary as each trajectory is about 1 GB
    if (!keepTracks)
      gSystem->Exec(
          TString::Format("rm -f %s/track%d.root", outputDir.data(), n));

    std::cout << "\n";
  }  // Loop over electrons

  delete bathtubField;

  fout->cd();

  startTree->Write();
  hZPos->Write();
  hRPos->Write();
  h2RZPos->Write();
  hPitchAngle->Write();
  hInitialAngle->Write();
  hZPosAcc->Write();
  hRPosAcc->Write();
  h2RZPosAcc->Write();
  hPitchAngleAcc->Write();
  hInitialAngleAcc->Write();

  fout->Close();
  delete fout;
}
