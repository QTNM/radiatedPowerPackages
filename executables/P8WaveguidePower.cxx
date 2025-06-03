/*
  P8WaveguidePower.cxx

  Check project 8 waveguide power calculations with Signal machinery
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
#include "TVector3.h"
#include "Waveguides/CircularWaveguide.h"
#include "Waveguides/RectangularWaveguide.h"

using namespace rad;

int main() {
  TString outputFile{"~/work/qtnm/outputs/P8WaveguidePower/plots.root"};
  TFile fout(outputFile, "recreate");

  // Electron kinematics
  const double eKE{30e3};                        // eV
  const double centralFreq{26e9};                // GHz
  const double eSpeed{GetSpeedFromKE(eKE, ME)};  // metres per second
  TVector3 eVel(eSpeed, 0, 0);
  // Determine appropriate magnetic field
  const double bMag{2 * M_PI * centralFreq *
                    (ME + eKE * TMath::Qe() / pow(TMath::C(), 2)) /
                    TMath::Qe()};
  std::cout << "Required magnetic field = " << bMag << " T\n";

  // Calculate radiated power
  const double radPower{CalcLarmorPower(eKE, bMag, M_PI_2)};  // Watts
  std::cout << "Radiated power = " << radPower * 1e15 << " fW\n";

  // Waveguide length in metres
  const double wgLength{10e-2};
  TVector3 probePos(0, 0, wgLength / 2);
  // WR42 dimensions
  const double wr42LongSide{10.668e-3};  // m
  const double wr42ShortSide{4.318e-3};  // m
  auto wgWR42 = new RectangularWaveguide(wr42LongSide, wr42ShortSide, wgLength);
  // Circular waveguide radius in metres
  const double circWgRadius{5e-3};
  auto wgCirc = new CircularWaveguide(circWgRadius, wgLength);

  // Define magnetic field
  auto field = new UniformField(bMag);

  TString trackFile{"~/work/qtnm/outputs/P8WaveguidePower/track.root"};

  const double simTime{0.5e-6};  // seconds

  // Electron gyroradius in metres
  const double gyroradius{GetGyroradius(
      eVel, field->evaluate_field_at_point(TVector3(0, 0, 0)), ME)};
  const double xMin{-4.5e-3};  // metres
  const double xMax{4.5e-3};   // metres
  const uint nXPnts{31};
  auto grPowerWR42 = new TGraph();
  setGraphAttr(grPowerWR42);
  grPowerWR42->SetTitle(
      "WR42 waveguide; X [mm]; Fraction of total power (both directions)");
  grPowerWR42->SetMarkerStyle(20);

  for (uint iX{0}; iX < nXPnts; iX++) {
    const double x{xMin + (xMax - xMin) * double(iX) / double(nXPnts - 1)};
    std::cout << "X = " << x * 1e3 << " mm\n";
    TVector3 ePos(x, -gyroradius, 0);
    const double simStepSize{1.0 / (centralFreq * 15)};  // seconds
    ElectronTrajectoryGen traj(trackFile, field, ePos, eVel, simStepSize,
                               simTime);
    const double sampleRate{1e9};                       // Hertz
    const double loFreq{centralFreq - sampleRate / 4};  // Hertz
    LocalOscillator lo(2 * M_PI * loFreq);
    Signal sig(trackFile, wgWR42, lo, sampleRate,
               Probe(probePos, WaveguideMode(1, 0, kTE)));
    auto grP{sig.GetVIPowerPeriodogram(1)};
    grP->SetTitle(Form("X = %.2f mm", x * 1e3));
    fout.cd();
    grP->Write(Form("grP%d", iX));

    double powerIntegral{0};  // Integral of power in watts
    for (int i{0}; i < grP->GetN(); i++) {
      powerIntegral += grP->GetPointY(i);
    }
    grPowerWR42->SetPoint(grPowerWR42->GetN(), x * 1e3,
                          powerIntegral * 2 / radPower);

    delete grP;
  }

  // Now do the same thing with the circular waveguide
  const double rhoMin{0};       // metres
  const double rhoMax{4.5e-3};  // metres
  const uint nRhoPnts{31};
  auto grPowerCirc = new TGraph();
  setGraphAttr(grPowerCirc);
  grPowerCirc->SetTitle(
      "Circular waveguide; #rho [mm]; Fraction of total power (both "
      "directions)");
  grPowerCirc->SetMarkerStyle(20);

  for (uint iRho{0}; iRho < nRhoPnts; iRho++) {
    const double rho{rhoMin +
                     (rhoMax - rhoMin) * double(iRho) / double(nRhoPnts - 1)};
    std::cout << "Rho = " << rho * 1e3 << " mm\n";
    TVector3 ePos(0, rho - gyroradius, 0);
    const double simStepSize{1.0 / (centralFreq * 15)};  // seconds
    ElectronTrajectoryGen traj(trackFile, field, ePos, eVel, simStepSize,
                               simTime);
    const double sampleRate{1e9};                       // Hertz
    const double loFreq{centralFreq - sampleRate / 4};  // Hertz
    LocalOscillator lo(2 * M_PI * loFreq);
    Signal sig(trackFile, wgCirc, lo, sampleRate,
               Probe(probePos, WaveguideMode(1, 1, kTE)));
    auto grP{sig.GetVIPowerPeriodogram(1)};
    grP->SetTitle(Form("#rho = %.2f mm", rho * 1e3));
    fout.cd();
    grP->Write(Form("grP%d", iRho));

    double powerIntegral{0};  // Integral of power in watts
    for (int i{0}; i < grP->GetN(); i++) {
      powerIntegral += grP->GetPointY(i);
    }
    grPowerCirc->SetPoint(grPowerCirc->GetN(), rho * 1e3,
                          powerIntegral * 2 / radPower);

    delete grP;
  }
  fout.cd();
  grPowerWR42->Write("grPowerWR42");
  grPowerCirc->Write("grPowerCirc");
  fout.Close();
  delete grPowerWR42;
  delete grPowerCirc;
  return 0;
}