/*
 * PlotRuddXSecs.cxx
 *
 * Executable for makign some verification plots of inelastics cross-sections
 * Created on: 04-02-2025
 */

#include <getopt.h>
#include <unistd.h>

#include <cmath>
#include <iostream>
#include <string>

#include "BasicFunctions/BasicFunctions.h"
#include "Scattering/ElasticScatter.h"
#include "Scattering/InelasticScatter.h"
#include "TAxis.h"
#include "TFile.h"
#include "TGraph.h"

using namespace rad;

int main(int argc, char* argv[]) {
  int opt{};
  std::string outputStr{" "};  // Location of output ROOT file

  while ((opt = getopt(argc, argv, ":o:")) != -1) {
    switch (opt) {
      case 'o':
        outputStr = optarg;
        break;
      case ':':
        std::cerr << "Option -" << static_cast<char>(optopt)
                  << " requires an argument.\n";
        return 1;
      case '?':
        std::cerr << "Unknown option: " << static_cast<char>(optopt) << "\n";
        return 1;
    }
  }

  std::cout << "Writing output to " << outputStr << "\n";

  // Create ROOT file
  TFile fT(outputStr.c_str(), "recreate");

  /////////////////////////////////////////////////////////////////////
  //////////////////// Ionisation Cross Sections //////////////////////
  /////////////////////////////////////////////////////////////////////
  // Create graphs for total ionisation cross-sections in centred units
  // This is what is plotted in the Rudd paper
  TGraph grHRudd;
  grHRudd.SetName("grHRudd");
  grHRudd.SetTitle("H");
  setGraphAttr(grHRudd);
  grHRudd.SetMarkerStyle(20);
  grHRudd.SetMarkerColor(kRed);
  grHRudd.SetLineColor(kRed);
  grHRudd.SetLineWidth(2);
  grHRudd.GetXaxis()->SetTitle("T - I [eV]");
  grHRudd.GetYaxis()->SetTitle("Cross Section [10^{-20} m^{2}]");
  TGraph grH2Rudd;
  grH2Rudd.SetName("grH2Rudd");
  grH2Rudd.SetTitle("H2");
  setGraphAttr(grH2Rudd);
  grH2Rudd.SetMarkerStyle(20);
  grH2Rudd.SetMarkerColor(kOrange + 1);
  grH2Rudd.SetLineColor(kOrange + 1);
  grH2Rudd.SetLineWidth(2);
  grH2Rudd.GetXaxis()->SetTitle("T - I [eV]");
  grH2Rudd.GetYaxis()->SetTitle("Cross Section [10^{-20} m^{2}]");
  TGraph grHeRudd;
  grHeRudd.SetName("grHeRudd");
  grHeRudd.SetTitle("He");
  setGraphAttr(grHeRudd);
  grHeRudd.SetMarkerStyle(20);
  grHeRudd.SetMarkerColor(kGreen + 2);
  grHeRudd.SetLineColor(kGreen + 2);
  grHeRudd.SetLineWidth(2);
  grHeRudd.GetXaxis()->SetTitle("T - I [eV]");
  grHeRudd.GetYaxis()->SetTitle("Cross Section [10^{-20} m^{2}]");

  // Now the non-centred plots which are more useful for our analysis
  TGraph grH;
  grH.SetName("grH");
  grH.SetTitle("H");
  setGraphAttr(grH);
  grH.SetMarkerStyle(20);
  grH.SetMarkerColor(kRed);
  grH.SetLineColor(kRed);
  grH.SetLineWidth(2);
  grH.GetXaxis()->SetTitle("T [eV]");
  grH.GetYaxis()->SetTitle("Cross Section [10^{-20} m^{2}]");
  TGraph grH2;
  grH2.SetName("grH2");
  grH2.SetTitle("H2");
  setGraphAttr(grH2);
  grH2.SetMarkerStyle(20);
  grH2.SetMarkerColor(kOrange + 1);
  grH2.SetLineColor(kOrange + 1);
  grH2.SetLineWidth(2);
  grH2.GetXaxis()->SetTitle("T - I [eV]");
  grH2.GetYaxis()->SetTitle("Cross Section [10^{-20} m^{2}]");
  TGraph grHe;
  grHe.SetName("grHe");
  grHe.SetTitle("He");
  setGraphAttr(grHe);
  grHe.SetMarkerStyle(20);
  grHe.SetMarkerColor(kGreen + 2);
  grHe.SetLineColor(kGreen + 2);
  grHe.SetLineWidth(2);
  grHe.GetXaxis()->SetTitle("T - I [eV]");
  grHe.GetYaxis()->SetTitle("Cross Section [10^{-20} m^{2}]");

  const double maxTMinusI{20000};  // Kinetic minus binding energy in eV
  const double minTMinusI{0.5};    // Kinetic minus binding energy in eV
  const unsigned int nPoints{400};
  double stepSize{log10(maxTMinusI / minTMinusI) / (nPoints - 1)};
  const double I_H{RYDBERG_EV};  // Binding energy of H in eV
  const double I_H2{15.43};      // Binding energy of H2 in eV
  const double I_He{24.59};      // Binding energy of He in eV
  const double maxT{20e3};
  const double minT{I_He + 1.0};
  double stepSizeT{log10(maxT / minT) / (nPoints - 1)};
  for (unsigned int iPnt{0}; iPnt < nPoints; ++iPnt) {
    double TMinusI{minTMinusI * pow(10, iPnt * stepSize)};
    InelasticScatter scatterH(TMinusI + I_H, Species::H);
    InelasticScatter scatterH2(TMinusI + I_H2, Species::H2);
    InelasticScatter scatterHe(TMinusI + I_He, Species::He);
    grHRudd.SetPoint(iPnt, TMinusI, scatterH.GetTotalXSec() * 1e20);
    grH2Rudd.SetPoint(iPnt, TMinusI, scatterH2.GetTotalXSec() * 1e20);
    grHeRudd.SetPoint(iPnt, TMinusI, scatterHe.GetTotalXSec() * 1e20);

    double T{minT * pow(10, iPnt * stepSizeT)};
    InelasticScatter scatterH_2(T, Species::H);
    InelasticScatter scatterH2_2(T, Species::H2);
    InelasticScatter scatterHe_2(T, Species::He);
    grH.SetPoint(iPnt, T, scatterH.GetTotalXSec() * 1e20);
    grH2.SetPoint(iPnt, T, scatterH2.GetTotalXSec() * 1e20);
    grHe.SetPoint(iPnt, T, scatterHe.GetTotalXSec() * 1e20);
  }

  // Plots of the differential cross section in terms of W
  TGraph grdW_H;
  setGraphAttr(grdW_H);
  grdW_H.SetName("grdW_H");
  grdW_H.SetTitle("H");
  grdW_H.GetXaxis()->SetTitle("W [eV]");
  grdW_H.GetYaxis()->SetTitle("d#sigma/dW [m^{2}/eV]");
  grdW_H.SetLineColor(kRed);
  TGraph grdW_H2;
  setGraphAttr(grdW_H2);
  grdW_H2.SetName("grdW_H2");
  grdW_H2.SetTitle("H2");
  grdW_H2.GetXaxis()->SetTitle("W [eV]");
  grdW_H2.GetYaxis()->SetTitle("d#sigma/dW [m^{2}/eV]");
  grdW_H2.SetLineColor(kOrange + 1);
  TGraph grdW_He;
  setGraphAttr(grdW_He);
  grdW_He.SetName("grdW_He");
  grdW_He.SetTitle("He");
  grdW_He.GetXaxis()->SetTitle("W [eV]");
  grdW_He.GetYaxis()->SetTitle("d#sigma/dW [m^{2}/eV]");
  grdW_He.SetLineColor(kGreen + 2);

  const double wMin{1e-1};
  const double wMax{1e2};
  const unsigned int nPntsW{200};
  const double stepSizeW{log10(wMax / wMin) / (nPntsW - 1)};
  const double chosenT{18.6e3};
  InelasticScatter scatterH(chosenT, Species::H);
  InelasticScatter scatterH2(chosenT, Species::H2);
  InelasticScatter scatterHe(chosenT, Species::He);
  for (size_t i{0}; i < nPntsW; i++) {
    const double W{wMin * pow(10, i * stepSizeW)};
    grdW_H.SetPoint(i, W, scatterH.GetSDCS_W(W));
    grdW_H2.SetPoint(i, W, scatterH2.GetSDCS_W(W));
    grdW_He.SetPoint(i, W, scatterHe.GetSDCS_W(W));
  }

  ///////////// Elastic cross sections ////////////////////
  // Make some plots of the total elastic cross-sections //
  /////////////////////////////////////////////////////////
  const double TMin{200};
  const double TMax{20e3};
  const unsigned int nPnts{200};
  TGraph grHElastic;
  setGraphAttr(grHElastic);
  grHElastic.SetName("grHElastic");
  grHElastic.SetTitle("H");
  grHElastic.GetXaxis()->SetTitle("T [eV]");
  grHElastic.GetYaxis()->SetTitle("Cross Section [10^{-20} m^{2}]");
  grHElastic.SetLineColor(kMagenta + 1);
  TGraph grHeElastic;
  setGraphAttr(grHeElastic);
  grHeElastic.SetName("grHeElastic");
  grHeElastic.SetTitle("He");
  grHeElastic.GetXaxis()->SetTitle("T [eV]");
  grHeElastic.GetYaxis()->SetTitle("Cross Section [10^{-20} m^{2}]");
  grHeElastic.SetLineColor(kCyan + 1);

  TGraph grHRutherford;
  setGraphAttr(grHRutherford);
  grHRutherford.SetName("grHRutherford");
  grHRutherford.SetTitle("H");
  grHRutherford.GetXaxis()->SetTitle("T [eV]");
  grHRutherford.GetYaxis()->SetTitle("Cross Section [10^{-20} m^{2}]");
  grHRutherford.SetLineColor(kMagenta + 1);
  grHRutherford.SetLineStyle(2);
  TGraph grHeRutherford;
  setGraphAttr(grHeRutherford);
  grHeRutherford.SetName("grHeRutherford");
  grHeRutherford.SetTitle("He");
  grHeRutherford.GetXaxis()->SetTitle("T [eV]");
  grHeRutherford.GetYaxis()->SetTitle("Cross Section [10^{-20} m^{2}]");
  grHeRutherford.SetLineColor(kCyan + 1);
  grHRutherford.SetLineStyle(2);

  for (size_t i{0}; i < nPnts; i++) {
    double T{TMin + i * (TMax - TMin) / (nPnts - 1)};
    ElasticScatter scatterH(T, 1, 1);
    ElasticScatter scatterHe(T, 2, 4);
    grHElastic.SetPoint(i, T, scatterH.GetTotalXSec() * 1e20);
    grHeElastic.SetPoint(i, T, scatterHe.GetTotalXSec() * 1e20);
    grHRutherford.SetPoint(i, T, scatterH.TotalRutherfordXSec() * 1e20);
    grHeRutherford.SetPoint(i, T, scatterHe.TotalRutherfordXSec() * 1e20);
  }

  // Plots of the momentum transfer as a function of the scattering angle
  TGraph grMomentumTransfer;
  setGraphAttr(grMomentumTransfer);
  const double EIncident{18e3};  // eV
  grMomentumTransfer.SetName("grMomentumTransfer");
  grMomentumTransfer.SetTitle("Momentum Transfer, T = 18 keV");
  grMomentumTransfer.GetXaxis()->SetTitle("#theta [degrees]");
  grMomentumTransfer.GetYaxis()->SetTitle("q [kg m/s]");
  grMomentumTransfer.SetLineColor(kBlue);
  grMomentumTransfer.SetLineWidth(3);
  TGraph grOutgoingE;
  setGraphAttr(grOutgoingE);
  grOutgoingE.SetName("grOutgoingE");
  grOutgoingE.SetTitle("Outgoing KE, T = 18 keV");
  grOutgoingE.GetXaxis()->SetTitle("#theta [degrees]");
  grOutgoingE.GetYaxis()->SetTitle("E [keV]");
  grOutgoingE.SetLineColor(kBlue);
  grOutgoingE.SetLineWidth(3);

  const double gamma{1 + EIncident / ME_EV};
  const double beta{sqrt(1 - 1 / pow(gamma, 2))};
  const double pIncident{gamma * beta * ME * TMath::C()};

  const double thetaMin{1e-5};
  const double thetaMax{179 * M_PI / 180};
  stepSize = log10(thetaMax / thetaMin) / (nPnts - 1);
  for (size_t i{0}; i < nPnts; i++) {
    double theta{thetaMin * pow(10, i * stepSize)};
    ElasticScatter scatter(EIncident, 1, 1);
    double q{2 * pIncident * sin(theta / 2)};
    grMomentumTransfer.SetPoint(i, theta * 180 / M_PI, q);
    grOutgoingE.SetPoint(i, theta * 180 / M_PI,
                         scatter.GetEnergyAfterScatter(theta) / 1e3);
  }

  // Generate a histogram of the angular scattering distribution
  TH1D hTheta{"hTheta", "Angular Scattering Distribution", 100, 0, 0.005};
  TH1D hE{"hE", "Outgoing energy", 100, EIncident - 100, EIncident};
  ElasticScatter scatter(EIncident, 1, 1);
  for (size_t i{0}; i < 1000; i++) {
    double theta{scatter.GetRandomScatteringAngle()};
    hTheta.Fill(theta);
    hE.Fill(scatter.GetEnergyAfterScatter(theta));
  }

  fT.cd();
  grH.Write();
  grH2.Write();
  grHe.Write();
  grHRudd.Write();
  grH2Rudd.Write();
  grHeRudd.Write();

  grdW_H.Write();
  grdW_H2.Write();
  grdW_He.Write();

  grHElastic.Write();
  grHeElastic.Write();
  grHRutherford.Write();
  grHeRutherford.Write();
  grMomentumTransfer.Write();
  grOutgoingE.Write();

  hTheta.Write();
  hE.Write();

  fT.Close();
  return 0;
}