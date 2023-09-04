/*
  Project8CCA.cxx

  Take a look at the mode frequencies that are present in the hypothesised
  Project 8 CCA cavity design
*/

#include <iostream>

#include "BasicFunctions/BasicFunctions.h"
#include "BasicFunctions/Constants.h"
#include "TFile.h"
#include "TGraph.h"
#include "TLatex.h"
#include "TMarker.h"
#include "TMultiGraph.h"
#include "TString.h"
#include "TText.h"
#include "Waveguides/CircularCavity.h"

using namespace rad;

int main(int argc, char* argv[]) {
  TString plotFile{"/home/sjones/work/qtnm/Cavities/Project8CCA/plots.root"};
  TFile fout(plotFile, "recreate");

  const double innerDiameter{1.42e-2};          // metres
  const double innerRadius{innerDiameter / 2};  // metres
  const double lengthMin{10e-2};                // metres
  const double lengthMax{15e-2};                // metres

  const double energyKr83mLine{32.1e3};  // eV
  const double bFieldMag{1.0};           // Tesla
  const double cycFreq{CalcCyclotronFreq(energyKr83mLine, bFieldMag)};
  std::cout << "Cyclotron frequency = " << cycFreq / 1e9 << " GHz\n";

  TVector3 probePosition(innerRadius, 0, 0);
  CircularCavity cavMin(innerRadius, lengthMin, probePosition);
  const double fTE011Min{cavMin.GetResonantModeF(CircularCavity::kTE, 0, 1, 1)};
  const double fTE012Min{cavMin.GetResonantModeF(CircularCavity::kTE, 0, 1, 2)};
  auto grCavMinTE01p = new TGraph();
  setGraphAttr(grCavMinTE01p);
  grCavMinTE01p->SetMarkerStyle(20);
  grCavMinTE01p->SetMarkerColor(kRed);
  grCavMinTE01p->SetPoint(grCavMinTE01p->GetN(), fTE011Min / 1e9, 0.5);
  grCavMinTE01p->SetPoint(grCavMinTE01p->GetN(), fTE012Min / 1e9, 0.5);
  auto tTE01pMin = new TText(fTE011Min / 1e9, 0.55, "TE01p");
  tTE01pMin->SetTextColor(kRed);
  fout.cd();
  tTE01pMin->Write("tTE01pMin");

  const double fTE111Min{cavMin.GetResonantModeF(CircularCavity::kTE, 1, 1, 1)};
  const double fTE112Min{cavMin.GetResonantModeF(CircularCavity::kTE, 1, 1, 2)};
  auto grCavMinTE11p = new TGraph();
  setGraphAttr(grCavMinTE11p);
  grCavMinTE11p->SetMarkerStyle(20);
  grCavMinTE11p->SetMarkerColor(kBlue);
  grCavMinTE11p->SetPoint(grCavMinTE11p->GetN(), fTE111Min / 1e9, 1);
  grCavMinTE11p->SetPoint(grCavMinTE11p->GetN(), fTE112Min / 1e9, 1);
  auto tTE11pMin = new TText(fTE111Min / 1e9, 1.05, "TE11p");
  tTE11pMin->SetTextColor(kBlue);
  fout.cd();
  tTE11pMin->Write("tTE11pMin");

  const double fTE211Min{cavMin.GetResonantModeF(CircularCavity::kTE, 2, 1, 1)};
  const double fTE212Min{cavMin.GetResonantModeF(CircularCavity::kTE, 2, 1, 2)};
  auto grCavMinTE21p = new TGraph();
  setGraphAttr(grCavMinTE21p);
  grCavMinTE21p->SetMarkerStyle(20);
  grCavMinTE21p->SetMarkerColor(kBlack);
  grCavMinTE21p->SetPoint(grCavMinTE21p->GetN(), fTE211Min / 1e9, 1);
  grCavMinTE21p->SetPoint(grCavMinTE21p->GetN(), fTE212Min / 1e9, 1);
  auto tTE21pMin = new TText(fTE211Min / 1e9, 1.05, "TE21p");
  tTE21pMin->SetTextColor(kBlack);
  fout.cd();
  tTE21pMin->Write("tTE21pMin");

  const double fTM010Min{cavMin.GetResonantModeF(CircularCavity::kTM, 0, 1, 0)};
  const double fTM011Min{cavMin.GetResonantModeF(CircularCavity::kTM, 0, 1, 1)};
  const double fTM012Min{cavMin.GetResonantModeF(CircularCavity::kTM, 0, 1, 2)};
  auto grCavMinTM01p = new TGraph();
  setGraphAttr(grCavMinTM01p);
  grCavMinTM01p->SetMarkerStyle(20);
  grCavMinTM01p->SetMarkerColor(kOrange + 1);
  grCavMinTM01p->SetPoint(grCavMinTM01p->GetN(), fTM010Min / 1e9, 1);
  grCavMinTM01p->SetPoint(grCavMinTM01p->GetN(), fTM011Min / 1e9, 1);
  grCavMinTM01p->SetPoint(grCavMinTM01p->GetN(), fTM012Min / 1e9, 1);
  auto tTM01pMin = new TText(fTM010Min / 1e9, 1.05, "TM01p");
  tTM01pMin->SetTextColor(kOrange + 1);
  fout.cd();
  tTM01pMin->Write("tTM01pMin");

  const double fTM110Min{cavMin.GetResonantModeF(CircularCavity::kTM, 1, 1, 0)};
  const double fTM111Min{cavMin.GetResonantModeF(CircularCavity::kTM, 1, 1, 1)};
  const double fTM112Min{cavMin.GetResonantModeF(CircularCavity::kTM, 1, 1, 2)};
  auto grCavMinTM11p = new TGraph();
  setGraphAttr(grCavMinTM11p);
  grCavMinTM11p->SetMarkerStyle(20);
  grCavMinTM11p->SetMarkerColor(kMagenta + 1);
  grCavMinTM11p->SetPoint(grCavMinTM11p->GetN(), fTM110Min / 1e9, 1);
  grCavMinTM11p->SetPoint(grCavMinTM11p->GetN(), fTM111Min / 1e9, 1);
  grCavMinTM11p->SetPoint(grCavMinTM11p->GetN(), fTM112Min / 1e9, 1);
  auto tTM11pMin = new TText(fTM110Min / 1e9, 1.05, "TM11p");
  tTM11pMin->SetTextColor(kMagenta + 1);
  fout.cd();
  tTM11pMin->Write("tTM11pMin");

  auto mgMin = new TMultiGraph("mgMin", "r = 7.2 mm, l = 100 mm; f [GHz]; ");
  mgMin->Add(grCavMinTE01p);
  mgMin->Add(grCavMinTE11p);
  mgMin->Add(grCavMinTE21p);
  mgMin->Add(grCavMinTM01p);
  mgMin->Add(grCavMinTM11p);

  std::cout << "TE011 = " << fTE011Min / 1e9
            << " GHz\tTE012 = " << fTE012Min / 1e9 << " GHz\n";
  std::cout << "TE111 = " << fTE111Min / 1e9
            << " GHz\tTE112 = " << fTE112Min / 1e9 << " GHz\n";
  std::cout << "TE211 = " << fTE211Min / 1e9
            << " GHz\tTE212 = " << fTE212Min / 1e9 << " GHz\n";
  std::cout << "TM010 = " << fTM010Min / 1e9
            << " GHz\tTM011 = " << fTM011Min / 1e9
            << " GHz\tTM012 = " << fTM012Min / 1e9 << " GHz\n";
  std::cout << "TM110 = " << fTM110Min / 1e9
            << " GHz\tTM111 = " << fTM111Min / 1e9
            << " GHz\tTM112 = " << fTM112Min / 1e9 << " GHz\n";
  fout.cd();
  mgMin->Write();

  std::cout << "\n";

  CircularCavity cavMax(innerRadius, lengthMax, probePosition);
  const double fTE011Max{cavMax.GetResonantModeF(CircularCavity::kTE, 0, 1, 1)};
  const double fTE012Max{cavMax.GetResonantModeF(CircularCavity::kTE, 0, 1, 2)};
  auto grCavMaxTE01p = new TGraph();
  setGraphAttr(grCavMaxTE01p);
  grCavMaxTE01p->SetMarkerStyle(20);
  grCavMaxTE01p->SetMarkerColor(kRed);
  grCavMaxTE01p->SetPoint(grCavMaxTE01p->GetN(), fTE011Max / 1e9, 0.5);
  grCavMaxTE01p->SetPoint(grCavMaxTE01p->GetN(), fTE012Max / 1e9, 0.5);
  auto tTE01pMax = new TText(fTE011Max / 1e9, 0.55, "TE01p");
  tTE01pMax->SetTextColor(kRed);
  fout.cd();
  tTE01pMax->Write("tTE01pMax");

  const double fTE111Max{cavMax.GetResonantModeF(CircularCavity::kTE, 1, 1, 1)};
  const double fTE112Max{cavMax.GetResonantModeF(CircularCavity::kTE, 1, 1, 2)};
  auto grCavMaxTE11p = new TGraph();
  setGraphAttr(grCavMaxTE11p);
  grCavMaxTE11p->SetMarkerStyle(20);
  grCavMaxTE11p->SetMarkerColor(kBlue);
  grCavMaxTE11p->SetPoint(grCavMaxTE11p->GetN(), fTE111Max / 1e9, 1);
  grCavMaxTE11p->SetPoint(grCavMaxTE11p->GetN(), fTE112Max / 1e9, 1);
  auto tTE11pMax = new TText(fTE111Max / 1e9, 1.05, "TE11p");
  tTE11pMax->SetTextColor(kBlue);
  fout.cd();
  tTE11pMax->Write("tTE11pMax");

  const double fTE211Max{cavMax.GetResonantModeF(CircularCavity::kTE, 2, 1, 1)};
  const double fTE212Max{cavMax.GetResonantModeF(CircularCavity::kTE, 2, 1, 2)};
  auto grCavMaxTE21p = new TGraph();
  setGraphAttr(grCavMaxTE21p);
  grCavMaxTE21p->SetMarkerStyle(20);
  grCavMaxTE21p->SetMarkerColor(kBlack);
  grCavMaxTE21p->SetPoint(grCavMaxTE21p->GetN(), fTE211Max / 1e9, 1);
  grCavMaxTE21p->SetPoint(grCavMaxTE21p->GetN(), fTE212Max / 1e9, 1);
  auto tTE21pMax = new TText(fTE211Max / 1e9, 1.05, "TE21p");
  tTE21pMax->SetTextColor(kBlack);
  fout.cd();
  tTE21pMax->Write("tTE21pMax");

  const double fTM010Max{cavMax.GetResonantModeF(CircularCavity::kTM, 0, 1, 0)};
  const double fTM011Max{cavMax.GetResonantModeF(CircularCavity::kTM, 0, 1, 1)};
  const double fTM012Max{cavMax.GetResonantModeF(CircularCavity::kTM, 0, 1, 2)};
  auto grCavMaxTM01p = new TGraph();
  setGraphAttr(grCavMaxTM01p);
  grCavMaxTM01p->SetMarkerStyle(20);
  grCavMaxTM01p->SetMarkerColor(kOrange + 1);
  grCavMaxTM01p->SetPoint(grCavMaxTM01p->GetN(), fTM010Max / 1e9, 1);
  grCavMaxTM01p->SetPoint(grCavMaxTM01p->GetN(), fTM011Max / 1e9, 1);
  grCavMaxTM01p->SetPoint(grCavMaxTM01p->GetN(), fTM012Max / 1e9, 1);
  auto tTM01pMax = new TText(fTM010Max / 1e9, 1.05, "TM01p");
  tTM01pMax->SetTextColor(kOrange + 1);
  fout.cd();
  tTM01pMax->Write("tTM01pMax");

  const double fTM110Max{cavMax.GetResonantModeF(CircularCavity::kTM, 1, 1, 0)};
  const double fTM111Max{cavMax.GetResonantModeF(CircularCavity::kTM, 1, 1, 1)};
  const double fTM112Max{cavMax.GetResonantModeF(CircularCavity::kTM, 1, 1, 2)};
  auto grCavMaxTM11p = new TGraph();
  setGraphAttr(grCavMaxTM11p);
  grCavMaxTM11p->SetMarkerStyle(20);
  grCavMaxTM11p->SetMarkerColor(kMagenta + 1);
  grCavMaxTM11p->SetPoint(grCavMaxTM11p->GetN(), fTM110Max / 1e9, 1);
  grCavMaxTM11p->SetPoint(grCavMaxTM11p->GetN(), fTM111Max / 1e9, 1);
  grCavMaxTM11p->SetPoint(grCavMaxTM11p->GetN(), fTM112Max / 1e9, 1);
  auto tTM11pMax = new TText(fTM110Max / 1e9, 1.05, "TM11p");
  tTM11pMax->SetTextColor(kMagenta + 1);
  fout.cd();
  tTM11pMax->Write("tTM11pMax");

  auto mgMax = new TMultiGraph("mgMax", "r = 7.2 mm, l = 150 mm; f [GHz]; ");
  mgMax->Add(grCavMaxTE01p);
  mgMax->Add(grCavMaxTE11p);
  mgMax->Add(grCavMaxTE21p);
  mgMax->Add(grCavMaxTM01p);
  mgMax->Add(grCavMaxTM11p);

  std::cout << "TE011 = " << fTE011Max / 1e9
            << " GHz\tTE012 = " << fTE012Max / 1e9 << " GHz\n";
  std::cout << "TE111 = " << fTE111Max / 1e9
            << " GHz\tTE112 = " << fTE112Max / 1e9 << " GHz\n";
  std::cout << "TE211 = " << fTE211Max / 1e9
            << " GHz\tTE212 = " << fTE212Max / 1e9 << " GHz\n";
  std::cout << "TM010 = " << fTM010Max / 1e9
            << " GHz\tTM011 = " << fTM011Max / 1e9
            << " GHz\tTM012 = " << fTM012Max / 1e9 << " GHz\n";
  std::cout << "TM110 = " << fTM110Max / 1e9
            << " GHz\tTM111 = " << fTM111Max / 1e9
            << " GHz\tTM112 = " << fTM112Max / 1e9 << " GHz\n";
  fout.cd();
  mgMax->Write();

  fout.Close();
  return 0;
}