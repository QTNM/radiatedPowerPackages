/*
  PhasedArrayVariableField.cxx

  Build phased arrays for a hypothetical large trap
*/

#include <getopt.h>

#include <cmath>
#include <iostream>

#include "Antennas/HalfWaveDipole.h"
#include "BasicFunctions/BasicFunctions.h"
#include "BasicFunctions/Constants.h"
#include "ElectronDynamics/QTNMFields.h"
#include "ElectronDynamics/TrajectoryGen.h"
#include "SignalProcessing/InducedVoltage.h"
#include "SignalProcessing/LocalOscillator.h"
#include "SignalProcessing/NoiseFunc.h"
#include "SignalProcessing/Signal.h"
#include "TFile.h"
#include "TGraph.h"
#include "TMath.h"
#include "TString.h"
#include "TTree.h"
#include "TVector3.h"

using namespace rad;

int main(int argc, char *argv[]) {
  int opt;

  std::string outputDir{" "};
  double bField{1.0};  // Tesla

  while ((opt = getopt(argc, argv, ":d:b:")) != -1) {
    switch (opt) {
      case 'd':
        outputDir = optarg;
        break;
      case 'b':
        bField = atof(optarg);
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
  if (bField <= 0) {
    std::cout << "Invalid B field value of " << bField << " T" << std::endl;
    exit(1);
  }

  std::cout << "Output directory is " << outputDir << std::endl;
  std::cout << "Selected B field is " << bField << " T" << std::endl;

  const clock_t begin_time = clock();

  TFile *fout = new TFile(Form("%s/output.root", outputDir.data()), "RECREATE");

  // Create a bathtub field
  const double vacuumFlangeDiameter{160e-3};
  const double trapDepth{5e-3};  // Tesla
  const double trapLength{1.0};  // metres
  const double iCoil{2 * trapDepth * (vacuumFlangeDiameter / 2.0) / MU0};
  BathtubField *field =
      new BathtubField(vacuumFlangeDiameter / 2.0, iCoil, -trapLength / 2.0,
                       trapLength / 2.0, TVector3(0, 0, bField));
  const double centralField{
      field->evaluate_field_at_point(TVector3(0, 0, 0)).Mag()};

  const double thetaBot{asin(
      sqrt(1.0 - trapDepth / (trapDepth + bField)))};  // Minimum trapping angle
  std::cout << "Minimum trapping pitch angle = "
            << thetaBot * 180.0 / TMath::Pi() << " degrees" << std::endl;

  // Electron kinematics
  const double electronKE{18600};  // eV
  const double electronSpeed{GetSpeedFromKE(electronKE, ME)};
  const double tau = 2 * R_E / (3 * TMath::C());
  TVector3 V0(electronSpeed, 0, 0);
  const double gyroradius{
      GetGyroradius(V0, field->evaluate_field_at_point(TVector3(0, 0, 0)), ME)};
  TVector3 X0(0, -gyroradius, 0);

  const double centralFreq{CalcCyclotronFreq(electronKE, centralField)};
  const double centralPeriod{1.0 / centralFreq};
  const double centralLambda{TMath::C() / centralFreq};
  std::cout << "Electron frequency, wavelength = " << centralFreq / 1e9
            << " GHz, " << centralLambda * 1e2 << " cm" << std::endl;

  // Setup the track writing
  const double simTime{55e-6};
  const double simStepSize{1e-12};

  TString trackFile{Form("%s/track.root", outputDir.data())};
  ElectronTrajectoryGen traj(trackFile, field, X0, V0, simStepSize, simTime,
                             0.0, tau);
  std::cout << "Generated electron track" << std::endl;

  // Generate the antenna array
  const double antennaZoneDiameter{122e-3};
  const double antennaZoneRadius{antennaZoneDiameter / 2.0};
  // Assuming a half wave dipole plus a little more for fitting these in
  const double antennaSize{(centralLambda / 2.0) + 1e-2};
  std::cout << "Antennas have a physical size of " << antennaSize * 1e3 << " mm"
            << std::endl;
  const int nAntennas{
      int(floor((TMath::Pi() * antennaZoneDiameter) / antennaSize))};
  std::cout << "We have " << nAntennas << " antennas" << std::endl;

  std::vector<IAntenna *> antennaArray;
  const double loadResistance{73.0};

  for (int iAnt = 0; iAnt < nAntennas; iAnt++) {
    double antennaAngle{2 * TMath::Pi() * double(iAnt) / double(nAntennas)};
    double timeShift{centralPeriod * (1.0 - double(iAnt) / double(nAntennas))};
    std::cout << "Dipole " << iAnt
              << ": Angle = " << antennaAngle * 180 / TMath::Pi()
              << " degrees. Required time shift is " << timeShift * 1e12
              << " ps" << std::endl;
    TVector3 antennaPoint(antennaZoneRadius * cos(antennaAngle),
                          antennaZoneRadius * sin(antennaAngle), 0.0);
    TVector3 antennaDirZ(-1 * sin(antennaAngle), cos(antennaAngle), 0.0);
    TVector3 antennaDirX(cos(antennaAngle), sin(antennaAngle), 0.0);
    HalfWaveDipole *antenna = new HalfWaveDipole(
        antennaPoint, antennaDirX, antennaDirZ, centralFreq, timeShift);
    antennaArray.push_back(antenna);

    InducedVoltage iv(trackFile, {antenna}, true);
    iv.GenerateVoltage(0, 5e-8);
    TGraph *grV{iv.GetVoltageGraph()};
    TGraph *grPower{iv.GetPowerPeriodogram(loadResistance)};

    fout->cd();
    grV->Write(Form("grV%d", iAnt));
    grPower->Write(Form("grPower%d", iAnt));

    delete grPower;
    delete grV;
    std::cout << "\n";
  }

  // Now create an InducedVoltage for the combined signal
  InducedVoltage ivCombined(trackFile, antennaArray, true);
  ivCombined.GenerateVoltage(0, 5e-8);
  TGraph *grVCombined{ivCombined.GetVoltageGraph()};
  TGraph *grPowerCombined{ivCombined.GetPowerPeriodogram(loadResistance)};

  fout->cd();
  grVCombined->Write("grVCombined");
  grPowerCombined->Write("grPowerCombined");
  delete grVCombined;
  delete grPowerCombined;

  // Now we have the short time stuff, looking at an actual processed signal
  const double tAcq{simTime - 1e-6};  // seconds
  const double sampleRate{750e6};     // Hertz
  const double loFreq{centralFreq - sampleRate / 4.0};
  LocalOscillator lo(2 * TMath::Pi() * loFreq);
  const double noiseTemp{4.0};  // Kelvin
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

  Signal sigNoNoise(trackFile, antennaArray, lo, sampleRate, {});
  std::cout << "Created the signal with no noise" << std::endl;
  TGraph *grVSigNoNoise = sigNoNoise.GetVITimeDomain();
  TGraph *grVSigNoNoisePgram = sigNoNoise.GetVIPowerPeriodogram(loadResistance);
  fout->cd();
  grVSigNoNoise->Write("grVSigNoNoise");
  grVSigNoNoisePgram->Write("grVSigNoNoisePgram");
  delete grVSigNoNoise;
  delete grVSigNoNoisePgram;

  const clock_t end_time = clock();
  std::cout << "Simulation time was "
            << float(end_time - begin_time) / CLOCKS_PER_SEC << " seconds"
            << std::endl;

  fout->Close();
  delete fout;
  return 0;
}
