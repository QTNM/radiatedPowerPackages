/*
 * testXSecs.cxx
 *
 * Executable for making verification plots of cross-sections
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
  grH2.GetXaxis()->SetTitle("T [eV]");
  grH2.GetYaxis()->SetTitle("Cross Section [10^{-20} m^{2}]");
  TGraph grHe;
  grHe.SetName("grHe");
  grHe.SetTitle("He");
  setGraphAttr(grHe);
  grHe.SetMarkerStyle(20);
  grHe.SetMarkerColor(kGreen + 2);
  grHe.SetLineColor(kGreen + 2);
  grHe.SetLineWidth(2);
  grHe.GetXaxis()->SetTitle("T [eV]");
  grHe.GetYaxis()->SetTitle("Cross Section [10^{-20} m^{2}]");

  const double maxTMinusI{20000};  // Kinetic minus binding energy in eV
  const double minTMinusI{0.5};    // Kinetic minus binding energy in eV
  const unsigned int nPoints{1000};
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
    grH.SetPoint(iPnt, T, scatterH_2.GetTotalXSec() * 1e20);
    grH2.SetPoint(iPnt, T, scatterH2_2.GetTotalXSec() * 1e20);
    grHe.SetPoint(iPnt, T, scatterHe_2.GetTotalXSec() * 1e20);
  }

  /////////////////////////////////////////////////////////////////////
  //////// Data from Shah et al. 1987 J. Phys. B 20, 3501   ///////////
  /////////////////////////////////////////////////////////////////////
  double energiesShah[106] = {
      14.6,   14.8,   15.0,   15.1,   15.2,   15.4,   15.6,   15.9,   16.1,
      16.4,   16.6,   16.9,   17.1,   17.4,   17.6,   17.9,   18.1,   18.4,
      18.7,   19.0,   19.3,   19.6,   20.0,   20.4,   20.9,   21.4,   22.0,
      22.6,   23.3,   24.0,   24.8,   25.6,   26.6,   27.3,   28.3,   29.3,
      30.5,   31.6,   32.8,   34.1,   35.4,   36.7,   38.1,   39.6,   41.2,
      42.9,   44.7,   46.6,   48.6,   50.7,   52.9,   55.2,   57.6,   60.1,
      63.0,   66.0,   69.0,   72.1,   75.5,   79.5,   84.0,   89.0,   94.0,
      102.0,  103.0,  113.0,  121.0,  130.2,  138.2,  148.2,  158.2,  168.2,
      178.2,  188.2,  198.2,  213.2,  228.2,  248.2,  268.2,  288.0,  317.9,
      347.9,  387.9,  427.9,  467.9,  508.2,  548.2,  598.2,  668.2,  748.2,
      818.2,  898.2,  998.2,  1100.0, 1200.0, 1300.0, 1506.7, 1662.7, 1848.1,
      1998.1, 2198.1, 2448.1, 2698.1, 2998.1, 3298.1, 3648.1};
  double xsecHShah[106] = {
      0.544e-21, 0.661e-21, 0.762e-21, 0.820e-21, 0.870e-21, 0.990e-21,
      1.08e-21,  1.25e-21,  1.37e-21,  1.45e-21,  1.63e-21,  1.68e-21,
      1.73e-21,  1.96e-21,  2.07e-21,  2.15e-21,  2.22e-21,  2.35e-21,
      2.50e-21,  2.61e-21,  2.75e-21,  2.81e-21,  2.93e-21,  3.11e-21,
      3.34e-21,  3.39e-21,  3.61e-21,  3.76e-21,  4.01e-21,  4.15e-21,
      4.30e-21,  4.44e-21,  4.57e-21,  4.75e-21,  4.95e-21,  5.01e-21,
      5.10e-21,  5.27e-21,  5.39e-21,  5.53e-21,  5.59e-21,  5.74e-21,
      5.83e-21,  5.89e-21,  6.02e-21,  6.07e-21,  6.08e-21,  6.23e-21,
      6.27e-21,  6.19e-21,  6.23e-21,  6.21e-21,  6.13e-21,  6.14e-21,
      6.11e-21,  6.11e-21,  6.01e-21,  5.96e-21,  5.91e-21,  5.84e-21,
      5.78e-21,  5.59e-21,  5.40e-21,  5.42e-21,  5.23e-21,  5.07e-21,
      5.05e-21,  4.83e-21,  4.62e-21,  4.55e-21,  4.43e-21,  4.28e-21,
      4.10e-21,  3.98e-21,  3.79e-21,  3.61e-21,  3.43e-21,  3.31e-21,
      3.03e-21,  2.84e-21,  2.66e-21,  2.50e-21,  2.31e-21,  2.15e-21,
      2.00e-21,  1.86e-21,  1.77e-21,  1.59e-21,  1.47e-21,  1.38e-21,
      1.26e-21,  1.13e-21,  1.05e-21,  0.982e-21, 0.914e-21, 0.807e-21,
      0.721e-21, 0.673e-21, 0.631e-21, 0.577e-21, 0.525e-21, 0.472e-21,
      0.437e-21, 0.403e-21, 0.370e-21, 0.339e-21};
  for (size_t i{0}; i < 106; i++) {
    xsecHShah[i] *= 1e20;
  }
  TGraph grHShah(106, energiesShah, xsecHShah);
  setGraphAttr(grHShah);
  grHShah.SetName("grHShah");
  grHShah.SetTitle("Shah (1987)");
  grHShah.GetXaxis()->SetTitle("T [eV]");
  grHShah.GetYaxis()->SetTitle("Cross Section [m^{2}]");
  grHShah.SetMarkerStyle(5);
  grHShah.SetMarkerColor(kYellow + 2);

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

  for (size_t i{0}; i < nPoints; i++) {
    double T{minT * pow(10, double(i) * stepSizeT)};
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

  grHShah.Write();

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