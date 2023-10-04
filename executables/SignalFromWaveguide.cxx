/*
  SignalFromWaveguide
*/

#include <cmath>
#include <iostream>

#include "BasicFunctions/BasicFunctions.h"
#include "BasicFunctions/Constants.h"
#include "ElectronDynamics/QTNMFields.h"
#include "ElectronDynamics/TrajectoryGen.h"
#include "SignalProcessing/LocalOscillator.h"
#include "SignalProcessing/Signal.h"
#include "TFile.h"
#include "TGraph.h"
#include "TString.h"
#include "TTreeReader.h"
#include "TTreeReaderValue.h"
#include "TVector3.h"
#include "Waveguides/CircularWaveguide.h"

using namespace rad;

int main() {
  TString plotFile{
      "/home/sjones/work/qtnm/outputs/SignalFromWaveguide/plots.root"};
  TFile fPlot(plotFile, "recreate");

  // Set up magnetic field and trap
  const double bBkg{0.7};                           // Tesla
  const double trapDepth{3.4e-3};                   // Tesla
  const double rLoop{0.02};                         // metres
  const double iLoop{2 * trapDepth * rLoop / MU0};  // Amps
  auto field = new HarmonicField(rLoop, iLoop, bBkg);
  TVector3 centralB{field->evaluate_field_at_point(TVector3(0, 0, 0))};

  // Electron kinematics
  const double eKE{18.6e3};
  const double eSpeed{GetSpeedFromKE(eKE, ME)};
  const double cycFreq{CalcCyclotronFreq(eKE, centralB.Mag())};
  std::cout << "Cyclotron frequency = " << cycFreq / 1e9 << " GHz\n";
  const double pitchAngleDeg{87.8};
  const double pitchAngleRad{pitchAngleDeg * M_PI / 180};
  TVector3 eVel(eSpeed * sin(pitchAngleRad), 0, eSpeed * cos(pitchAngleRad));
  const double rg{GetGyroradius(eVel, centralB, ME)};
  TVector3 ePos(0, rg, 0);

  TString trackFile{
      "/home/sjones/work/qtnm/outputs/SignalFromWaveguide/track.root"};
  const double simTime{2e-6};                    // seconds
  const double simStepSize{1 / (cycFreq * 15)};  // seconds
  std::cout << "Generating electron trajectory\n";
  ElectronTrajectoryGen traj(trackFile, field, ePos, eVel, simStepSize, simTime,
                             true);

  // Set up a couple of waveguides with different probe positions
  const double wgRadius{7.14e-3};  // metres
  const double wgLength{0.08};     // metres
  TVector3 probePos1(wgRadius * 0.5, 0, 0);
  auto wg1 = new CircularWaveguide(wgRadius, wgLength, probePos1);
  TVector3 probePos2(0, wgRadius * 0.5, 0);
  auto wg2 = new CircularWaveguide(wgRadius, wgLength, probePos2);

  // Set up signals for these two things
  const double sampleRate{1e9};
  const double loFreq{cycFreq - 220e6};
  std::cout << "Local oscillator frequency = " << loFreq / 1e9 << " GHz\n";
  LocalOscillator lo(loFreq * TMath::TwoPi());
  Signal sig1(trackFile, wg1, lo, sampleRate);
  std::cout << "Finished generating first signal\n";
  Signal sig2(trackFile, wg2, lo, sampleRate);
  std::cout << "Finished generating second signal\n";

  auto grV1T{sig1.GetVITimeDomain()};
  grV1T->SetLineColor(kBlack);
  auto grV1P{sig1.GetVIPowerPeriodogram(1)};
  grV1P->SetLineColor(kBlack);

  auto grV2T{sig2.GetVITimeDomain()};
  grV2T->SetLineColor(kRed);
  auto grV2P{sig2.GetVIPowerPeriodogram(1)};
  grV2P->SetLineColor(kRed);

  fPlot.cd();
  grV1T->Write("grV1T");
  grV1P->Write("grV1P");
  grV2T->Write("grV2T");
  grV2P->Write("grV2P");

  fPlot.Close();
  return 0;
}