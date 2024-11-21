/*
  TestWgTheory.cxx
*/

#include <getopt.h>
#include <unistd.h>

#include <boost/lexical_cast.hpp>
#include <boost/uuid/uuid.hpp>
#include <boost/uuid/uuid_generators.hpp>
#include <boost/uuid/uuid_io.hpp>
#include <cmath>
#include <ctime>
#include <filesystem>
#include <iostream>
#include <random>
#include <string>

#include "BasicFunctions/BasicFunctions.h"
#include "BasicFunctions/Constants.h"
#include "ElectronDynamics/QTNMFields.h"
#include "ElectronDynamics/TrajectoryGen.h"
#include "H5Cpp.h"
#include "Scattering/ElasticScatter.h"
#include "Scattering/InelasticScatter.h"
#include "SignalProcessing/LocalOscillator.h"
#include "SignalProcessing/Signal.h"
#include "TFile.h"
#include "TString.h"
#include "TSystem.h"
#include "TTree.h"
#include "TVector3.h"
#include "Waveguides/CircularWaveguide.h"
#include "Waveguides/Probe.h"
#include "Waveguides/RectangularWaveguide.h"
#include "Waveguides/WaveguideMode.h"

using namespace rad;
using std::cout;
using std::endl;

/// @brief Function for generating a random file string
/// @return A unique string which can be used for file names
std::string make_uuid() {
  return boost::lexical_cast<std::string>((boost::uuids::random_generator())());
}

int main(int argc, char* argv[]) {
  int opt{};
  std::string outputStemStr{" "};  // Directory in which to store files
  double maxSimTime{10e-6};        // seconds

  while ((opt = getopt(argc, argv, ":o:t:")) != -1) {
    switch (opt) {
      case 'o':
        outputStemStr = optarg;
        break;
      case 't':
        maxSimTime = boost::lexical_cast<double>(optarg);
        break;
      case ':':
        cout << "Option -" << static_cast<char>(optopt)
             << " requires an argument." << endl;
        return 1;
      case '?':
        cout << "Unknown option: " << static_cast<char>(optopt) << endl;
        return 1;
    }
  }

  cout << "Writing output files to " << outputStemStr << "\n";
  cout << "Simulation time = " << maxSimTime * 1e6 << " us\n";

  // Check if the chosen output directory exists
  bool dirExists{std::filesystem::is_directory(outputStemStr)};
  // If it doesn't exist then attempt to create the directory
  if (!dirExists) {
    bool dirCreated{std::filesystem::create_directories(outputStemStr)};
    if (!dirCreated) {
      cout << "Failed to create the output directory. Exiting.\n";
      exit(1);
    }
  }

  // Create a directory for the temporary track files
  std::string trackDir{outputStemStr + "/tracks"};
  bool trackDirExists{std::filesystem::is_directory(trackDir)};
  if (!trackDirExists) {
    bool trackDirCreated{std::filesystem::create_directories(trackDir)};
    if (!trackDirCreated) {
      cout << "Failed to create the track directory. Exiting.\n";
      exit(1);
    }
  }

  TString outputDir{outputStemStr};

  const double frequency{26e9};  // 26 GHz
  const double energy{18.6e3};   // 30 keV
  const double gamma{energy / ME_EV + 1};
  const double reqBField{2 * M_PI * frequency * ME * gamma / TMath::Qe()};
  cout << "Required B field = " << reqBField << " Tesla\n";
  const double freeSpacePower{CalcLarmorPower(energy, reqBField, M_PI / 2)};
  cout << "Free space power = " << freeSpacePower * 1e15 << " fW\n";

  // Define uniform background field
  auto field{new UniformField(reqBField)};
  const unsigned int nRadialSteps{21};

  // Define the waveguide
  const double wr42Width{10.7e-3};
  const double wr42Height{4.3e-3};

  auto wg{new RectangularWaveguide(wr42Width, wr42Height, 10e-2)};
  const double eSpeed{GetSpeedFromKE(energy, ME)};
  const double gyroradius{
      GetGyroradius(TVector3(eSpeed, 0, 0),
                    field->evaluate_field_at_point(TVector3(0, 0, 0)), ME)};

  auto grPower{new TGraph()};
  auto grPowerSignal{new TGraph()};

  WaveguideMode modeTE10(1, 0, kTE);
  Probe probeTE10(TVector3(0, 0, 5e-2), WaveguideMode(1, 0, kTE));

  TFile plots("plots.root", "RECREATE");

  // Generate a series of electrons
  for (size_t iStep{0}; iStep < nRadialSteps; iStep++) {
    const double x{-4.5e-3 + double(iStep) / double(nRadialSteps - 1) * 9e-3};
    TVector3 initPos{x + gyroradius, 0, 0};
    TVector3 initVel{0, eSpeed, 0};
    cout << "Generating electron at x = " << x * 1e3 << " mm\n";
    // Now generate the trajectory
    std::string trackFileExt{make_uuid()};
    TString trackFile{trackDir + Form("/track_%s.root", trackFileExt.data())};
    const double simStepSize{1e-12};  // seconds
    ElectronTrajectoryGen traj(trackFile, field, initPos, initVel, simStepSize,
                               1e-6);

    const double sampleRate{1e9};                     // Hz
    const double loFreq{frequency - sampleRate / 4};  // Hz
    // Define local oscillator
    LocalOscillator lo(2 * M_PI * loFreq);
    Signal sigPlus(trackFile, wg, lo, sampleRate, probeTE10);
    auto grVPlus{sigPlus.GetVITimeDomain()};

    double avgPower{0};
    const double modeImp{wg->GetModeImpedance(modeTE10, 2 * M_PI * frequency)};
    const double k_c{wg->GetCutoffWavenumber(modeTE10)};
    const double omega{2 * M_PI * frequency};
    const double beta{sqrt(pow(omega / TMath::C(), 2) - k_c * k_c)};
    const double modePn{M_PI * M_PI * beta * MU0 * omega *
                        (wr42Height * wr42Height) /
                        (2 * wr42Height * wr42Width * pow(k_c, 4))};
    cout << "Mode impedance = " << modeImp << endl;
    cout << "Pn = " << modePn << endl;

    const double ANormCalc{2 * k_c * k_c / (M_PI * MU0 * omega) *
                           sqrt(wr42Height * wr42Width / pow(wr42Height, 2))};
    const double integralPlus{
        wg->GetEFieldIntegral(modeTE10, 2 * M_PI * frequency, 1.0, 100, true)};
    const double normPlus{1 / sqrt(2 * integralPlus)};
    cout << "Norm+ = " << normPlus << endl;

    // Open the track file
    TFile fin(trackFile, "READ");
    TTree* trackTree = (TTree*)fin.Get("tree");
    double time{};
    double xPos{}, yPos{}, zPos{};
    double xVel{}, yVel{}, zVel{};
    trackTree->SetBranchAddress("time", &time);
    trackTree->SetBranchAddress("xPos", &xPos);
    trackTree->SetBranchAddress("yPos", &yPos);
    trackTree->SetBranchAddress("zPos", &zPos);
    trackTree->SetBranchAddress("xVel", &xVel);
    trackTree->SetBranchAddress("yVel", &yVel);
    trackTree->SetBranchAddress("zVel", &zVel);

    auto grAmpPozar{new TGraph()};
    grAmpPozar->SetTitle(
        Form("X = %.2f mm; Time [s]; Field amplitude", x * 1e3));
    auto grAmpJackson{new TGraph()};
    grAmpJackson->SetTitle(
        Form("X = %.2f mm; Time [s]; Field amplitude", x * 1e3));
    auto grPowerPozar{new TGraph()};
    grPowerPozar->SetTitle(Form("X = %.2f mm; Time [s]; P [W]", x * 1e3));
    for (int ii{0}; ii < trackTree->GetEntries(); ii++) {
      trackTree->GetEntry(ii);
      TVector3 thePos{xPos, yPos, zPos};
      TVector3 theVel{xVel, yVel, zVel};
      TVector3 modeField{wg->GetModeEField(thePos, modeTE10, normPlus,
                                           2 * M_PI * frequency, true)};

      TVector3 J{-TMath::Qe() * theVel};
      const double ATE10{-1 / modePn * modeField.Dot(J)};
      const double ampJackson{
          -modeImp / 2 *
          J.Dot(wg->GetModeEField(thePos, modeTE10, ANormCalc,
                                  2 * M_PI * frequency, true))};
      const double powerIntegralTE10{
          M_PI * M_PI * beta * MU0 * omega * (wr42Height * wr42Height) /
          (4 * wr42Height * wr42Width * pow(k_c, 4)) * pow(ATE10, 2)};
      grAmpPozar->SetPoint(ii, time, ATE10);
      grPowerPozar->SetPoint(ii, time, powerIntegralTE10);
      grAmpJackson->SetPoint(ii, time, ampJackson);
      avgPower += powerIntegralTE10;
    }

    double avgPowerSignal{0};
    for (int i{0}; i < grVPlus->GetN(); i++) {
      avgPowerSignal += grVPlus->GetPointY(i) * grVPlus->GetPointY(i);
    }
    avgPowerSignal /= grVPlus->GetN();
    avgPowerSignal /= modeImp;
    avgPowerSignal *= 2;
    grPowerSignal->SetPoint(iStep, x * 1e3, avgPowerSignal / freeSpacePower);

    avgPower /= trackTree->GetEntries();
    avgPower *= 2;
    cout << "Average power = " << avgPower * 1e15 << " fW" << endl;
    grPower->SetPoint(iStep, x * 1e3, avgPower / freeSpacePower);

    fin.Close();
    // Delete the track file
    std::filesystem::remove(trackFile.Data());

    plots.cd();
    grAmpPozar->Write(Form("ampPozar_%d", int(iStep)));
    grAmpJackson->Write(Form("ampJackson_%d", int(iStep)));
    grPowerPozar->Write(Form("powerPozar_%d", int(iStep)));
  }

  grPower->SetMarkerStyle(20);
  grPower->SetTitle("WR42; X [mm]; Fraction of free space power");
  grPower->Write("power");
  grPowerSignal->SetMarkerStyle(20);
  grPowerSignal->SetTitle("WR42; X [mm]; Fraction of free space power");
  grPowerSignal->Write("powerSignal");

  //////////////////////////////////////////////////////////////////////////////
  /// Now look at the case with the circular waveguide /////////////////////////
  //////////////////////////////////////////////////////////////////////////////
  const double wgRadius{5e-3};   // meters
  const double wgLength{10e-2};  // meters

  auto wgCirc{new CircularWaveguide(wgRadius, wgLength)};
  WaveguideMode modeTE11(1, 1, kTE);
  const double circModeImp{
      wgCirc->GetModeImpedance(modeTE11, 2 * M_PI * frequency)};
  cout << "Impedance of TE11 mode = " << circModeImp << endl;

  auto grPowerSigCirc_P{new TGraph()};
  grPowerSigCirc_P->SetTitle("TE11+; X [mm]; Fraction of free space power");
  setGraphAttr(grPowerSigCirc_P);
  auto grPowerSigCirc_M{new TGraph()};
  grPowerSigCirc_M->SetTitle("TE11-; X [mm]; Fraction of free space power");
  setGraphAttr(grPowerSigCirc_M);
  auto grPowerSigCirc_TM01{new TGraph()};
  grPowerSigCirc_TM01->SetTitle("TM01; X [mm]; Fraction of free space power");
  setGraphAttr(grPowerSigCirc_TM01);
  // Generate a series of electrons
  for (size_t iStep{0}; iStep < 15; iStep++) {
    const double x{0 + double(iStep) / double(15 - 1) * 4.5e-3};
    TVector3 initPos{x + gyroradius, 0, 0};
    TVector3 initVel{0, eSpeed, 0};
    cout << "Generating electron at x = " << x * 1e3 << " mm\n";
    // Now generate the trajectory
    std::string trackFileExt{make_uuid()};
    TString trackFile{trackDir + Form("/track_%s.root", trackFileExt.data())};
    const double simStepSize{1e-12};  // seconds
    ElectronTrajectoryGen traj(trackFile, field, initPos, initVel, simStepSize,
                               5e-7);

    const double sampleRate{1e9};                     // Hz
    const double loFreq{frequency - sampleRate / 4};  // Hz
    // Define local oscillator
    LocalOscillator lo(2 * M_PI * loFreq);

    // Get the signals
    // TE11 mode readouts
    // First polarisation
    Signal sigTE11_P_Z1(
        trackFile, wgCirc, lo, sampleRate,
        Probe(TVector3(0, 0, 5e-2), WaveguideMode(1, 1, kTE), true));
    auto grTE11_P_Z1{sigTE11_P_Z1.GetVITimeDomain()};
    Signal sigTE11_P_Z2(
        trackFile, wgCirc, lo, sampleRate,
        Probe(TVector3(0, 0, -5e-2), WaveguideMode(1, 1, kTE), true));
    auto grTE11_P_Z2{sigTE11_P_Z2.GetVITimeDomain()};
    plots.cd();
    grTE11_P_Z1->Write(Form("TE11_P_Z1_%d", int(iStep)));
    grTE11_P_Z2->Write(Form("TE11_P_Z2_%d", int(iStep)));
    // Second polarisation
    Signal sigTE11_M_Z1(
        trackFile, wgCirc, lo, sampleRate,
        Probe(TVector3(0, 0, 5e-2), WaveguideMode(1, 1, kTE), false));
    auto grTE11_M_Z1{sigTE11_M_Z1.GetVITimeDomain()};
    Signal sigTE11_M_Z2(
        trackFile, wgCirc, lo, sampleRate,
        Probe(TVector3(0, 0, -5e-2), WaveguideMode(1, 1, kTE), false));
    auto grTE11_M_Z2{sigTE11_M_Z2.GetVITimeDomain()};
    plots.cd();
    grTE11_M_Z1->Write(Form("TE11_M_Z1_%d", int(iStep)));
    grTE11_M_Z2->Write(Form("TE11_M_Z2_%d", int(iStep)));

    // TM01 readouts
    Signal sigTM01_Z1(trackFile, wgCirc, lo, sampleRate,
                      Probe(TVector3(0, 0, 5e-2), WaveguideMode(0, 1, kTM)));
    auto grVTM01_Z1{sigTM01_Z1.GetVITimeDomain()};
    Signal sigTM01_Z2(trackFile, wgCirc, lo, sampleRate,
                      Probe(TVector3(0, 0, -5e-2), WaveguideMode(0, 1, kTM)));
    auto grVTM01_Z2{sigTM01_Z2.GetVITimeDomain()};
    plots.cd();
    grVTM01_Z1->Write(Form("TM01_Z1_%d", int(iStep)));
    grVTM01_Z2->Write(Form("TM01_Z2_%d", int(iStep)));

    double modeImpTM01{wgCirc->GetModeImpedance(WaveguideMode(0, 1, kTM),
                                                2 * M_PI * frequency)};

    double avgPowerTM01{0};
    for (int i{0}; i < grVTM01_Z1->GetN(); i++) {
      avgPowerTM01 += grVTM01_Z1->GetPointY(i) * grVTM01_Z1->GetPointY(i);
      avgPowerTM01 += grVTM01_Z2->GetPointY(i) * grVTM01_Z2->GetPointY(i);
    }
    avgPowerTM01 /= (grVTM01_Z1->GetN() * modeImpTM01);
    grPowerSigCirc_TM01->SetPoint(iStep, x * 1e3,
                                  avgPowerTM01 / freeSpacePower);

    double avgPowerTE11_P{0};
    for (int i{0}; i < grTE11_P_Z1->GetN(); i++) {
      avgPowerTE11_P += pow(grTE11_P_Z1->GetPointY(i), 2);
      avgPowerTE11_P += pow(grTE11_P_Z2->GetPointY(i), 2);
    }
    avgPowerTE11_P /= (grTE11_P_Z1->GetN() * circModeImp);
    grPowerSigCirc_P->SetPoint(iStep, x * 1e3, avgPowerTE11_P / freeSpacePower);

    double avgPowerTE11_M{0};
    for (int i{0}; i < grTE11_M_Z1->GetN(); i++) {
      avgPowerTE11_M += pow(grTE11_M_Z1->GetPointY(i), 2);
      avgPowerTE11_M += pow(grTE11_M_Z2->GetPointY(i), 2);
    }
    avgPowerTE11_M /= (grTE11_M_Z1->GetN() * circModeImp);
    grPowerSigCirc_M->SetPoint(iStep, x * 1e3, avgPowerTE11_M / freeSpacePower);

    // Delete the track file
    std::filesystem::remove(trackFile.Data());
  }
  plots.cd();
  grPowerSigCirc_P->Write("powerCirc_TE11_P");
  grPowerSigCirc_M->Write("powerCirc_TE11_M");
  grPowerSigCirc_TM01->Write("powerCirc_TM01");

  plots.Close();

  return 0;
}