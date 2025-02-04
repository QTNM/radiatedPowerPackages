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

  // Create graphs for total cross-sections
  TGraph grH2;
  grH2.SetName("grH2");
  grH2.SetTitle("H2");
  grH2.SetMarkerStyle(20);
  grH2.SetMarkerColor(kOrange + 1);
  grH2.SetLineColor(kOrange + 1);
  grH2.SetLineWidth(3);
  grH2.GetXaxis()->SetTitle("T - I [eV]");
  grH2.GetYaxis()->SetTitle("Cross Section [10^{-20} m^{2}]");

  TGraph grHe;
  grHe.SetName("grHe");
  grHe.SetTitle("He");
  grHe.SetMarkerStyle(20);
  grHe.SetMarkerColor(kGreen + 2);
  grHe.SetLineColor(kGreen + 2);
  grHe.SetLineWidth(3);
  grHe.GetXaxis()->SetTitle("T - I [eV]");
  grHe.GetYaxis()->SetTitle("Cross Section [10^{-20} m^{2}]");

  const double maxTMinusI{20000};  // Kinetic minus binding energy in eV
  const double minTMinusI{0.5};    // Kinetic minus binding energy in eV
  const unsigned int nPoints{200};
  const double stepSize{log10(maxTMinusI / minTMinusI) / (nPoints - 1)};
  const double I_H2{15.43};  // Binding energy of H2 in eV
  const double I_He{24.59};  // Binding energy of He in eV
  for (unsigned int iPnt{0}; iPnt < nPoints; ++iPnt) {
    double tMinusI{minTMinusI * pow(10, iPnt * stepSize)};
    InelasticScatter scatterH2(tMinusI + I_H2, Species::H2);
    InelasticScatter scatterHe(tMinusI + I_He, Species::He);
    grH2.SetPoint(iPnt, tMinusI, scatterH2.GetTotalXSec() * 1e20);
    grHe.SetPoint(iPnt, tMinusI, scatterHe.GetTotalXSec() * 1e20);
  }

  fT.cd();
  grH2.Write();
  grHe.Write();

  fT.Close();
  return 0;
}