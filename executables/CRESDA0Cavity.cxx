/*
  CRESDA0Cavity.cxx

  Make mode plot for hypothetical CRESDA0 cavity
*/

#include <iostream>

#include "BasicFunctions/BasicFunctions.h"
#include "TFile.h"
#include "TGraph.h"
#include "TLine.h"
#include "TMultiGraph.h"
#include "TString.h"
#include "TText.h"
#include "Waveguides/CircularCavity.h"

using namespace rad;

int main() {
  TString plotFile{"/home/sjones/work/qtnm/Cavities/CRESDA0Cavity/plots.root"};
  TFile fout(plotFile, "recreate");

  // Define cavity
  const double tElec{18.6e3};                                 // eV
  const double bFieldMag{0.7};                                // T
  const double cavityRadius{5e-3};                            // metres
  const double cycFreq{CalcCyclotronFreq(tElec, bFieldMag)};  // Hz
  const double cycFreq0{CalcCyclotronFreq(0, bFieldMag)};     // Hz
  std::cout << "Cyclotron frequency = " << cycFreq / 1e9 << " GHz\n";
  std::cout << "Cyclotron frequency (0 eV) = " << cycFreq0 / 1e9 << " GHz\n";
  const double p11Prime{GetBesselPrimeZero(1, 1)};
  const double d5mm_TE111{1 * TMath::Pi() /
                          sqrt(pow(TMath::TwoPi() * cycFreq / TMath::C(), 2) -
                               pow(p11Prime / cavityRadius, 2))};
  TVector3 probePosition(cavityRadius, 0, 0);
  CircularCavity cav(cavityRadius, d5mm_TE111, probePosition);

  // Calculate resonant mode frequencies
  ////////////// TE //////////////
  const double fTE011{cav.GetResonantModeF(CircularCavity::kTE, 0, 1, 1)};
  const double fTE012{cav.GetResonantModeF(CircularCavity::kTE, 0, 1, 2)};
  auto grCavTE01p = new TGraph();
  setGraphAttr(grCavTE01p);
  grCavTE01p->SetMarkerStyle(20);
  grCavTE01p->SetMarkerColor(kRed);
  grCavTE01p->SetPoint(grCavTE01p->GetN(), fTE011 / 1e9, 1);
  grCavTE01p->SetPoint(grCavTE01p->GetN(), fTE012 / 1e9, 1);
  auto tTE01p = new TText(fTE011 / 1e9, 0.95, "TE01p");
  tTE01p->SetTextColor(kRed);
  fout.cd();
  tTE01p->Write("tTE01p");

  const double fTE111{cav.GetResonantModeF(CircularCavity::kTE, 1, 1, 1)};
  const double fTE112{cav.GetResonantModeF(CircularCavity::kTE, 1, 1, 2)};
  auto grCavTE11p = new TGraph();
  setGraphAttr(grCavTE11p);
  grCavTE11p->SetMarkerStyle(20);
  grCavTE11p->SetMarkerColor(kBlue);
  grCavTE11p->SetPoint(grCavTE11p->GetN(), fTE111 / 1e9, 1);
  grCavTE11p->SetPoint(grCavTE11p->GetN(), fTE112 / 1e9, 1);
  auto tTE11p = new TText(fTE111 / 1e9, 0.95, "TE11p");
  tTE11p->SetTextColor(kBlue);
  fout.cd();
  tTE11p->Write("tTE11p");

  const double fTE211{cav.GetResonantModeF(CircularCavity::kTE, 2, 1, 1)};
  const double fTE212{cav.GetResonantModeF(CircularCavity::kTE, 2, 1, 2)};
  auto grCavTE21p = new TGraph();
  setGraphAttr(grCavTE21p);
  grCavTE21p->SetMarkerStyle(20);
  grCavTE21p->SetMarkerColor(kBlack);
  grCavTE21p->SetPoint(grCavTE21p->GetN(), fTE211 / 1e9, 1);
  grCavTE21p->SetPoint(grCavTE21p->GetN(), fTE212 / 1e9, 1);
  auto tTE21p = new TText(fTE211 / 1e9, 0.95, "TE21p");
  tTE21p->SetTextColor(kBlack);
  fout.cd();
  tTE21p->Write("tTE21p");

  ////////////// TM //////////////
  const double fTM010{cav.GetResonantModeF(CircularCavity::kTM, 0, 1, 0)};
  const double fTM011{cav.GetResonantModeF(CircularCavity::kTM, 0, 1, 1)};
  const double fTM012{cav.GetResonantModeF(CircularCavity::kTM, 0, 1, 2)};
  auto grCavTM01p = new TGraph();
  setGraphAttr(grCavTM01p);
  grCavTM01p->SetMarkerStyle(20);
  grCavTM01p->SetMarkerColor(kOrange + 1);
  grCavTM01p->SetPoint(grCavTM01p->GetN(), fTM010 / 1e9, 0.5);
  grCavTM01p->SetPoint(grCavTM01p->GetN(), fTM011 / 1e9, 0.5);
  grCavTM01p->SetPoint(grCavTM01p->GetN(), fTM012 / 1e9, 0.5);
  auto tTM01p = new TText(fTM010 / 1e9, 0.55, "TM01p");
  tTM01p->SetTextColor(kOrange + 1);
  fout.cd();
  tTM01p->Write("tTM01p");

  const double fTM110{cav.GetResonantModeF(CircularCavity::kTM, 1, 1, 0)};
  const double fTM111{cav.GetResonantModeF(CircularCavity::kTM, 1, 1, 1)};
  const double fTM112{cav.GetResonantModeF(CircularCavity::kTM, 1, 1, 2)};
  auto grCavTM11p = new TGraph();
  setGraphAttr(grCavTM11p);
  grCavTM11p->SetMarkerStyle(20);
  grCavTM11p->SetMarkerColor(kMagenta + 1);
  grCavTM11p->SetPoint(grCavTM11p->GetN(), fTM110 / 1e9, 0.5);
  grCavTM11p->SetPoint(grCavTM11p->GetN(), fTM111 / 1e9, 0.5);
  grCavTM11p->SetPoint(grCavTM11p->GetN(), fTM112 / 1e9, 0.5);
  auto tTM11p = new TText(fTM110 / 1e9, 0.55, "TM11p");
  tTM11p->SetTextColor(kMagenta + 1);
  fout.cd();
  tTM11p->Write("tTM11p");

  const double fTM210{cav.GetResonantModeF(CircularCavity::kTM, 2, 1, 0)};
  const double fTM211{cav.GetResonantModeF(CircularCavity::kTM, 2, 1, 1)};
  const double fTM212{cav.GetResonantModeF(CircularCavity::kTM, 2, 1, 2)};
  auto grCavTM21p = new TGraph();
  setGraphAttr(grCavTM21p);
  grCavTM21p->SetMarkerStyle(20);
  grCavTM21p->SetMarkerColor(kGreen + 1);
  grCavTM21p->SetPoint(grCavTM21p->GetN(), fTM210 / 1e9, 0.5);
  grCavTM21p->SetPoint(grCavTM21p->GetN(), fTM211 / 1e9, 0.5);
  grCavTM21p->SetPoint(grCavTM21p->GetN(), fTM212 / 1e9, 0.5);
  auto tTM21p = new TText(fTM210 / 1e9 - 1, 0.55, "TM21p");
  tTM21p->SetTextColor(kGreen + 1);
  fout.cd();
  tTM21p->Write("tTM21p");

  auto mg = new TMultiGraph(
      "mg", Form("r = 5 mm, l = %.1f mm; f [GHz]; ", d5mm_TE111 * 1e3));
  mg->Add(grCavTE01p);
  mg->Add(grCavTE11p);
  mg->Add(grCavTE21p);
  mg->Add(grCavTM01p);
  mg->Add(grCavTM11p);
  mg->Add(grCavTM21p);

  TLine lMin(cycFreq / 1e9, 0.45, cycFreq / 1e9, 1);
  lMin.SetLineWidth(2);
  lMin.SetLineStyle(2);
  lMin.SetLineColor(kCyan + 1);

  TLine lMax(cycFreq0 / 1e9, 0.45, cycFreq0 / 1e9, 1);
  lMax.SetLineWidth(2);
  lMax.SetLineStyle(2);
  lMax.SetLineColor(kCyan + 1);

  fout.cd();
  mg->Write();
  lMin.Write("lMin");
  lMax.Write("lMax");
  fout.Close();

  return 0;
}