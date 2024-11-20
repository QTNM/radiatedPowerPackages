/*
  P8HarmonicTrap.cxx

  Look at recreating the harmonic trap analysis from arXiv:1901.02844

  S. Jones 05/09/23
*/

#include <getopt.h>
#include <unistd.h>

#include <cmath>
#include <iostream>

#include "BasicFunctions/BasicFunctions.h"
#include "BasicFunctions/Constants.h"
#include "ElectronDynamics/QTNMFields.h"
#include "ElectronDynamics/TrajectoryGen.h"
#include "SignalProcessing/LocalOscillator.h"
#include "SignalProcessing/Signal.h"
#include "TF1.h"
#include "TFile.h"
#include "TH2.h"
#include "TMath.h"
#include "TString.h"
#include "TTree.h"
#include "Waveguides/CircularWaveguide.h"

using namespace rad;

double HarmonicFieldFunc(double* x, double* par) {
  double z{x[0] / 100};
  double b0{par[0]};
  double l0{par[1]};
  return b0 * (1 + z * z / (l0 * l0));
}

double FindMainPeakPower(TGraph* gr, double searchFreq, double freqTolerance) {
  const double pntsPerFreq{double(gr->GetN()) / gr->GetPointX(gr->GetN() - 1)};
  const int startPnt{int(round(pntsPerFreq * (searchFreq - freqTolerance)))};
  const int endPnt{int(round(pntsPerFreq * (searchFreq + freqTolerance)))};

  // First find the highest point
  int peakInd{0};
  double peakHeight{-DBL_MAX};
  for (int i{startPnt}; i <= endPnt; i++) {
    if (gr->GetPointY(i) > peakHeight) {
      peakHeight = gr->GetPointY(i);
      peakInd = i;
    }
  }

  return gr->GetPointY(peakInd - 2) + gr->GetPointY(peakInd - 1) +
         gr->GetPointY(peakInd) + gr->GetPointY(peakInd + 1) +
         gr->GetPointY(peakInd + 2);
}

int main(int argc, char* argv[]) {
  TString plotFile{"/home/sjones/work/qtnm/outputs/P8HarmonicTrap/plots.root"};
  TFile fPlots(plotFile, "recreate");

  // Define the field
  const double l0{20e-2};        // m
  const double rCoil{1.2e-2};    // m
  const double trapDepth{4e-3};  // Tesla
  const double iCoils{2 * rCoil * trapDepth / MU0};
  const double bkgField{1.004};
  const double b0{bkgField - trapDepth};  // Tesla
  auto field = new HarmonicField(rCoil, iCoils, bkgField);
  TVector3 centralField{field->evaluate_field_at_point(TVector3(0, 0, 0))};

  const uint nZPnts{400};
  auto grProf = new TGraph();
  setGraphAttr(grProf);
  grProf->SetTitle("; z [cm]; |B| [T]");
  for (uint i{0}; i < nZPnts; i++) {
    const double z{-0.07 + 0.14 * double(i) / double(nZPnts - 1)};
    grProf->SetPoint(grProf->GetN(), z * 100,
                     -field->evaluate_field_at_point(TVector3(0, 0, z)).Z());
  }

  const double zFuncLimit{sqrt(trapDepth) * l0};
  auto fShape = new TF1("fShape", HarmonicFieldFunc, -zFuncLimit * 100,
                        zFuncLimit * 100, 2);
  fShape->SetParameter(0, b0);
  fShape->SetParameter(1, l0);

  // Electron kinematics (the fixed parts anyway)
  const double eKE{30e3};                        // eV
  const double eSpeed{GetSpeedFromKE(eKE, ME)};  // m/s
  const double beta{eSpeed / TMath::C()};
  std::cout << "Beta = " << beta << std::endl;
  const double gamma{1 / sqrt(1 - beta * beta)};
  const double rg{GetGyroradius(TVector3(eSpeed, 0, 0), centralField, ME)};
  TVector3 x0(0, rg, 0);
  std::cout << "x0 = " << x0.Y() * 1e3 << " mm\n";

  // Define the waveguide
  const double wgRadius{5e-3};            // metres
  const double wgLength{8e-2};            // metres
  TVector3 probePos(0, 0, wgLength / 2);  // Probe at the end of the guide
  auto wg = new CircularWaveguide(wgRadius, wgLength);
  const double cutoffTE11{wg->GetCutoffFrequency(WaveguideMode(1, 1, kTE))};
  const double cutoffTM01{wg->GetCutoffFrequency(WaveguideMode(0, 1, kTM))};
  std::cout << "Cutoff TE11 = " << cutoffTE11 / 1e9 << " GHz\n";
  std::cout << "Cutoff TM01 = " << cutoffTM01 / 1e9 << " GHz\n";

  const uint nzMaxPnts{41};
  const double zMaxMax{0.6e-2};  // m
  const double zMaxMin{0};       // m
  auto grMainPeak = new TGraph();
  setGraphAttr(grMainPeak);
  grMainPeak->SetTitle("Main peak; z_{max} [cm]; Relative power");
  grMainPeak->SetLineWidth(3);
  grMainPeak->SetLineColor(kBlue);
  grMainPeak->SetMarkerColor(kBlue);
  auto gr1stOrder = new TGraph();
  setGraphAttr(gr1stOrder);
  gr1stOrder->SetTitle("1st order sideband; z_{max} [cm]; Relative power");
  gr1stOrder->SetLineWidth(3);
  gr1stOrder->SetLineColor(kRed);
  gr1stOrder->SetMarkerColor(kRed);
  auto gr2ndOrder = new TGraph();
  setGraphAttr(gr2ndOrder);
  gr2ndOrder->SetTitle("2nd order sideband; z_{max} [cm]; Relative power");
  gr2ndOrder->SetLineWidth(3);
  gr2ndOrder->SetLineColor(kGreen + 1);
  gr2ndOrder->SetMarkerColor(kGreen + 1);

  for (uint iZ{0}; iZ < nzMaxPnts; iZ++) {
    const double zMax{zMaxMin +
                      (zMaxMax - zMaxMin) * double(iZ) / double(nzMaxPnts - 1)};
    const double thetaBotRad{atan2(l0, zMax)};
    const double thetaBotDeg{thetaBotRad * 180 / M_PI};
    const double avgCycFreq{(TMath::Qe() * b0 / (gamma * ME * 2 * M_PI)) *
                            (1 + zMax * zMax / (2 * l0 * l0))};
    const double axFreq{eSpeed * sin(thetaBotRad) / (l0 * 2 * M_PI)};
    std::cout << "z_max = " << zMax * 100 << " cm\ttheta_bot = " << thetaBotDeg
              << " degrees\tAverage freq. = " << avgCycFreq / 1e9
              << " GHz\tAxial freq. = " << axFreq / 1e6 << " MHz\n";

    // Simulate trajectory
    TVector3 v0(eSpeed * sin(thetaBotRad), 0, eSpeed * cos(thetaBotRad));
    const double simTime{1e-6};
    const double simTimeStep{1 / (15 * 27e9)};
    TString trackFile{
        "/home/sjones/work/qtnm/outputs/P8HarmonicTrap/track.root"};
    ElectronTrajectoryGen traj(trackFile, field, x0, v0, simTimeStep, simTime,
                               false);
    const double sampleRate{1.5e9};
    // Set up the local oscillator
    const double loFreq{avgCycFreq - sampleRate / 4};
    LocalOscillator lo(2 * M_PI * loFreq);
    Signal sig(trackFile, wg, lo, sampleRate,
               Probe(probePos, WaveguideMode(1, 1, kTE)));
    auto grV{sig.GetVIPowerPeriodogram(1)};
    grV->SetTitle(Form("z_{max} = %.1f cm", zMax * 100));

    fPlots.cd();
    grV->Write(Form("grV_%d", iZ));

    double peakPower{FindMainPeakPower(grV, sampleRate / 4, 5e6)};
    double order1Power1{FindMainPeakPower(grV, sampleRate / 4 + 90e6, 15e6)};
    double order1Power2{FindMainPeakPower(grV, sampleRate / 4 - 90e6, 15e6)};
    double order2Power1{
        FindMainPeakPower(grV, sampleRate / 4 + 2 * 90e6, 20e6)};
    double order2Power2{
        FindMainPeakPower(grV, sampleRate / 4 - 2 * 90e6, 20e6)};
    grMainPeak->SetPoint(grMainPeak->GetN(), zMax * 100, peakPower);
    gr1stOrder->SetPoint(gr1stOrder->GetN(), zMax * 100,
                         (order1Power1 + order1Power2) / 2);
    gr2ndOrder->SetPoint(gr2ndOrder->GetN(), zMax * 100,
                         (order2Power1 + order2Power2) / 2);
  }
  fPlots.cd();
  grProf->Write("grProf");
  fShape->Write();
  grMainPeak->Write("grMainPeak");
  gr1stOrder->Write("gr1stOrder");
  gr2ndOrder->Write("gr2ndOrder");

  fPlots.Close();
  return 0;
}