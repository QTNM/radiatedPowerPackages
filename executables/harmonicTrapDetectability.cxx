/// harmonicTrapDetectability.cxx

#include <getopt.h>

#include <cmath>
#include <iostream>

#include "Antennas/HalfWaveDipole.h"
#include "Antennas/IAntenna.h"
#include "BasicFunctions/Constants.h"
#include "ElectronDynamics/BorisSolver.h"
#include "ElectronDynamics/QTNMFields.h"
#include "ElectronDynamics/TrajectoryGen.h"
#include "SignalProcessing/InducedVoltage.h"
#include "SignalProcessing/LocalOscillator.h"
#include "SignalProcessing/NoiseFunc.h"
#include "SignalProcessing/Signal.h"
#include "TFile.h"
#include "TGraph.h"
#include "TMath.h"
#include "TMultiGraph.h"
#include "TRandom3.h"
#include "TSpline.h"
#include "TString.h"
#include "TTree.h"
#include "TVector3.h"

using namespace rad;

int main(int argc, char* argv[]) {
  int opt;
  std::string outputDir = " ";
  double pitchAngle = 90.0;
  double radialOffset = 0.0;
  int nDipoles = 33;

  while ((opt = getopt(argc, argv, ":d:p:r:n:")) != -1) {
    switch (opt) {
      case 'd':
        outputDir = optarg;
        break;
      case 'p':
        pitchAngle = atof(optarg);
        break;
      case 'r':
        radialOffset = atof(optarg);
        break;
      case 'n':
        nDipoles = atoi(optarg);
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
  if (outputDir == " ") {
    std::cout << "Must specify output file with -d" << std::endl;
    exit(1);
  }
  if (nDipoles <= 0) {
    std::cout << "Invalid number of dipoles (" << nDipoles
              << ") entered with -n" << std::endl;
  }

  std::cout << "Output directory is " << outputDir << std::endl;
  std::cout << "Pitch angle is " << pitchAngle << " degrees" << std::endl;
  std::cout << "Radial offset is " << (radialOffset * 1000) << " mm"
            << std::endl;
  std::cout << "Number of dipoles is " << nDipoles << std::endl;

  TString outputFile = Form("%s/output_r%.3f_p%.2f_n%d.root", outputDir.data(),
                            radialOffset, pitchAngle, nDipoles);
  std::cout << "Output file is " << outputFile << std::endl;

  TFile* fout = new TFile(outputFile, "RECREATE");

  const double coilRadius = 0.01;                           // metres
  const double bkgField = 1.0;                              // Tesla
  const double trapDepth = 4e-3;                            // Tesla
  const double ICoil = 2.0 * trapDepth * coilRadius / MU0;  // Amps

  HarmonicField* field = new HarmonicField(coilRadius, ICoil, bkgField);
  double centralField = field->evaluate_field_at_point(TVector3(0, 0, 0)).Mag();

  TGraph* grBMag = new TGraph();
  setGraphAttr(grBMag);
  grBMag->SetTitle("Harmonic trap; z [m]; |B| [T]");
  for (int i = 0; i < 321; i++) {
    double thisZ = -0.1 + 0.2 * double(i) / double(321);
    TVector3 BField = field->evaluate_field_at_point(TVector3(0, 0, thisZ));
    grBMag->SetPoint(grBMag->GetN(), thisZ, BField.Mag());
  }
  fout->cd();
  grBMag->Write("grBMag");

  // Start out with a 90 degree endpoint electron
  const double pitchAngleRad = pitchAngle * TMath::Pi() / 180.0;
  const double TElec = 18600;  // eV
  const double gamma = TElec * TMath::Qe() / (ME * TMath::C() * TMath::C()) + 1;
  const double betaSq = 1 - 1 / pow(gamma, 2);
  const double V0 = sqrt(betaSq) * TMath::C();
  const double tau = 2 * R_E / (3 * TMath::C());

  const double centralFrequency =
      (1.0 / (2 * TMath::Pi())) * TMath::Qe() * centralField / (gamma * ME);
  const double centralPeriod = 1.0 / centralFrequency;
  std::cout << "Central frequency = " << centralFrequency << " Hz" << std::endl;

  TVector3 vInitial(V0 * sin(pitchAngleRad), 0, V0 * cos(pitchAngleRad));
  const double gyroradius =
      gamma * ME * vInitial.X() / (TMath::Qe() * centralField);
  TVector3 X0(0, gyroradius + radialOffset, 0);
  std::cout << "Gyroradius = " << (gyroradius * 1000) << " mm" << std::endl;

  // Simulation time
  const double simTime = 55e-6;
  const double simStepSize = 1.5e-12;
  int nTimeSteps = simTime / simStepSize;

  // Specification for the antennas
  const double antennaRadius = 0.03;
  const double antennaLowerBandwidth = 26e9;
  const double antennaUpperBandwidth = 28e9;
  const double centralWavelength = TMath::C() / centralFrequency;
  const double boreCircumference = 2 * TMath::Pi() * antennaRadius;

  TString trackFilePath =
      Form("%s/track_r%.3f_p%.2f_n%d.root", outputDir.data(), radialOffset,
           pitchAngle, nDipoles);
  ElectronTrajectoryGen traj(trackFilePath, field, X0, vInitial, simStepSize,
                             simTime, 0.0, tau);

  // Firstly create our dipoles and add them to the array
  std::vector<IAntenna*> antennaArray;
  for (int iDip = 0; iDip < nDipoles; iDip++) {
    double antennaAngle1 = 2 * TMath::Pi() * double(iDip) / double(nDipoles);
    double reqShift = centralPeriod * double(iDip) / double(nDipoles);
    std::cout << "Dipole " << iDip
              << ": Angle = " << (antennaAngle1 * 180 / TMath::Pi())
              << " degrees. Required time shift is " << reqShift * 1e12 << " ps"
              << std::endl;

    // Antenna specifications
    TVector3 antennaPoint1(antennaRadius * TMath::Cos(antennaAngle1),
                           antennaRadius * TMath::Sin(antennaAngle1), 0.0);
    TVector3 antennaDirZ1(-1 * TMath::Sin(antennaAngle1),
                          TMath::Cos(antennaAngle1), 0.0);
    TVector3 antennaDirX1(TMath::Cos(antennaAngle1), TMath::Sin(antennaAngle1),
                          0.0);
    HalfWaveDipole* antenna = new HalfWaveDipole(
        antennaPoint1, antennaDirX1, antennaDirZ1, centralFrequency, reqShift);
    antenna->SetBandwidth(antennaLowerBandwidth, antennaUpperBandwidth);
    antennaArray.push_back(antenna);
  }

  // Now produce the voltage
  // Define the InducedVoltage as containing all the antennas in a single array

  // All sorts of signal processing stuff
  const double noiseTemp = 4.0;        // Kelvin
  const double loadResistance = 73.0;  // Ohms
  const double tAcq = simTime - 1e-6;  // seconds
  const double sampleRate = 750e6;     // Hz
  const double dmFreq = 26.7e9;        // Hz
  LocalOscillator lo(2 * TMath::Pi() * dmFreq);
  GaussianNoise noise1(noiseTemp, loadResistance);

  InducedVoltage iv(trackFilePath, antennaArray, true);
  Signal sig(trackFilePath, antennaArray, lo, sampleRate, {noise1});
  Signal sigNoNoise(trackFilePath, antennaArray, lo, sampleRate, {});

  TGraph* grVSig = sig.GetVITimeDomain();
  TGraph* grVSigPgram = sig.GetVIPowerPeriodogram(loadResistance);
  TGraph* grVSigNoNoise = sigNoNoise.GetVITimeDomain();
  TGraph* grVSigNoNoisePgram = sigNoNoise.GetVIPowerPeriodogram(loadResistance);

  fout->cd();
  grVSig->Write("grVSig");
  grVSigPgram->Write("grVSigPgram");
  grVSigNoNoise->Write("grVSigNoNoise");
  grVSigNoNoisePgram->Write("grVSigNoNoisePgram");

  delete grVSig;
  delete grVSigNoNoise;
  delete grVSigPgram;
  delete grVSigNoNoisePgram;

  fout->Close();
  delete fout;
  return 0;
}
