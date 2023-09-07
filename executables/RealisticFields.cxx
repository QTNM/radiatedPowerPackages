/// RealisticFields.cxx

#include <getopt.h>

#include <iostream>
#include <string>

#include "Antennas/HalfWaveDipole.h"
#include "Antennas/IAntenna.h"
#include "BasicFunctions/BasicFunctions.h"
#include "BasicFunctions/Constants.h"
#include "ElectronDynamics/ComsolFields.h"
#include "ElectronDynamics/QTNMFields.h"
#include "ElectronDynamics/TrajectoryGen.h"
#include "SignalProcessing/InducedVoltage.h"
#include "SignalProcessing/LocalOscillator.h"
#include "SignalProcessing/NoiseFunc.h"
#include "SignalProcessing/Signal.h"
#include "TFile.h"
#include "TGraph.h"
#include "TString.h"
#include "TTree.h"
#include "TVector3.h"

using namespace rad;

int main(int argc, char *argv[]) {
  int opt;

  std::string outputDir{" "};
  double desiredField{0.94};  // Tesla
  int nAntennas{1};
  unsigned int magnetSetup{1};
  double pitchAngle{90.0};
  double noiseTemp{0};

  while ((opt = getopt(argc, argv, ":d:b:n:m:p:t:")) != -1) {
    switch (opt) {
      case 'd':
        outputDir = optarg;
        break;
      case 'b':
        desiredField = atof(optarg);
        break;
      case 'n':
        nAntennas = atoi(optarg);
        break;
      case 'm':
        magnetSetup = atoi(optarg);
        break;
      case 'p':
        pitchAngle = atof(optarg);
        break;
      case 't':
        noiseTemp = atof(optarg);
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

  if (desiredField <= 0) {
    std::cout << "Invalid B field value of " << desiredField << " T\n";
    exit(1);
  }

  if (nAntennas <= 0) {
    std::cout << "Invalid number of antennas entered\n";
  }

  std::cout << "Using " << nAntennas << " antennas.\n";
  std::cout << "Chose a B field of " << desiredField << " T.\n";

  const double trapDepth{5e-3};  // Tesla
  // Choose magnet setup
  std::string fieldFileStem{
      "/home/sjones/work/qtnm/radiatedPowerPackages/files/"};
  std::string filePath{};
  double coilRadius{0};
  double antennaRadius{0};
  if (magnetSetup == 1) {
    // Option 1 is 110mm bore
    filePath = fieldFileStem + "fieldmap_50A_110mmx800mm.csv";
    coilRadius = 50e-3;
    antennaRadius = 28.5e-3;
    std::cout << "We are running the 110mm magnet bore.\n";
  } else if (magnetSetup == 2) {
    // Or we have the 170mm bore
    filePath = fieldFileStem + "fieldmap_50A_170mmx800mm.csv";
    coilRadius = 80e-3;
    antennaRadius = 55e-3;
    std::cout << "We are running the 170mm magnet bore.\n";
  }

  const double coilCurrent{2 * trapDepth * coilRadius / MU0};

  // Setup about fields
  double scaleFactor{desiredField / 0.94};
  ComsolHarmonicField *harm =
      new ComsolHarmonicField(coilRadius, coilCurrent, filePath, scaleFactor);

  const clock_t begin_time = clock();

  // Electron kinematics
  const double centralField{
      harm->evaluate_field_at_point(TVector3(0, 0.001, 0)).Mag()};
  const double electronKE{18600};  // eV
  const double electronSpeed{GetSpeedFromKE(electronKE, ME)};
  const double tau = 2 * R_E / (3 * TMath::C());
  TVector3 v0(electronSpeed, 0, 0);
  const double gyroradius{GetGyroradius(
      v0, harm->evaluate_field_at_point(TVector3(0, 0.001, 0)), ME)};
  TVector3 x0(0, -gyroradius, 0);

  TFile *fout = new TFile(
      Form("%s/output_%.2f.root", outputDir.data(), pitchAngle), "RECREATE");

  const double centralFreq{CalcCyclotronFreq(electronKE, centralField)};
  const double centralPeriod{1.0 / centralFreq};
  const double centralLambda{TMath::C() / centralFreq};
  std::cout << "Electron frequency, wavelength = " << centralFreq / 1e9
            << " GHz, " << centralLambda * 1e2 << " cm" << std::endl;

  // Setup the track writing
  const double simTime{55e-6};
  const double simStepSize{1e-12};
  TString trackFile{Form("%s/track_%.2f.root", outputDir.data(), pitchAngle)};
  ElectronTrajectoryGen traj(trackFile, harm, x0, v0, simStepSize, simTime);
  std::cout << "Generated trajectory\n";

  const double loadResistance{73.0};
  std::vector<IAntenna *> antennaArray;
  for (int iAnt = 0; iAnt < nAntennas; iAnt++) {
    double antennaAngle{2 * TMath::Pi() * double(iAnt) / double(nAntennas)};
    double timeShift{centralPeriod * (1.0 - double(iAnt) / double(nAntennas))};
    std::cout << "Dipole " << iAnt
              << ": Angle = " << antennaAngle * 180 / TMath::Pi()
              << " degrees. Required time shift is " << timeShift * 1e12
              << " ps" << std::endl;
    TVector3 antennaPoint(antennaRadius * cos(antennaAngle),
                          antennaRadius * sin(antennaAngle), 0.0);
    TVector3 antennaDirZ(-1 * sin(antennaAngle), cos(antennaAngle), 0.0);
    TVector3 antennaDirX(cos(antennaAngle), sin(antennaAngle), 0.0);
    HalfWaveDipole *antenna = new HalfWaveDipole(
        antennaPoint, antennaDirX, antennaDirZ, centralFreq, timeShift);
    antennaArray.push_back(antenna);
  }

  // Now we have the short time stuff, looking at an actual processed signal
  const double tAcq{simTime - 1e-6};  // seconds
  const double sampleRate{750e6};     // Hertz
  const double loFreq{centralFreq - sampleRate / 4.0};
  LocalOscillator lo(2 * TMath::Pi() * loFreq);
  GaussianNoise noiseFunc(noiseTemp, loadResistance);

  Signal sig(trackFile, antennaArray, lo, sampleRate, {noiseFunc});
  std::cout << "Created the signal with noise" << std::endl;
  TGraph *grVSig = sig.GetVITimeDomain();
  TGraph *grVSigPgram = sig.GetVIPowerPeriodogram(loadResistance);
  fout->cd();
  grVSig->Write("grVSig");
  grVSigPgram->Write("grVSigPgram");
  delete grVSig;
  delete grVSigPgram;

  const clock_t end_time = clock();
  std::cout << "Simulation time was "
            << float(end_time - begin_time) / CLOCKS_PER_SEC << " seconds"
            << std::endl;

  fout->Close();
  delete fout;
  return 0;
}
