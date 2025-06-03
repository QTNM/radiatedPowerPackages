/*
  LargeVolumeCavityP8.cxx

  Make plots of the modes for the large volume Project 8 cavity
  Details can be found in arXiv:2311.16415
*/

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
  TString plotFile{
      "/home/sjones/work/qtnm/Cavities/LargeVolumeCavityP8/plots.root"};
  TFile fout(plotFile, "recreate");

  const double cavDiameter{1.1};            // metres
  const double cavRadius{cavDiameter / 2};  // metres
  const double cavHeight{11};               // metres

  TVector3 probePosition(cavRadius, 0, 0);
  CircularCavity cavity(cavRadius, cavHeight, probePosition);

  // TE01p modes
  const double fTE011{cavity.GetResonantModeF(CircularCavity::kTE, 0, 1, 1)};
  const double fTE012{cavity.GetResonantModeF(CircularCavity::kTE, 0, 1, 2)};
  const double fTE013{cavity.GetResonantModeF(CircularCavity::kTE, 0, 1, 3)};
  const double fTE014{cavity.GetResonantModeF(CircularCavity::kTE, 0, 1, 4)};
  auto grTE01p = new TGraph();
  setGraphAttr(grTE01p);
  grTE01p->SetMarkerStyle(20);
  grTE01p->SetMarkerColor(kRed);
  grTE01p->SetPoint(grTE01p->GetN(), fTE011 / 1e6, 0.5);
  grTE01p->SetPoint(grTE01p->GetN(), fTE012 / 1e6, 0.5);
  grTE01p->SetPoint(grTE01p->GetN(), fTE013 / 1e6, 0.5);
  grTE01p->SetPoint(grTE01p->GetN(), fTE014 / 1e6, 0.5);
  auto tTE01p = new TText(fTE011 / 1e6, 0.55, "TE01p");
  tTE01p->SetTextColor(kRed);
  fout.cd();
  tTE01p->Write("tTE01p");

  // TE11p modes
  const double fTE111{cavity.GetResonantModeF(CircularCavity::kTE, 1, 1, 1)};
  const double fTE112{cavity.GetResonantModeF(CircularCavity::kTE, 1, 1, 2)};
  const double fTE113{cavity.GetResonantModeF(CircularCavity::kTE, 1, 1, 3)};
  const double fTE114{cavity.GetResonantModeF(CircularCavity::kTE, 1, 1, 4)};
  auto grTE11p = new TGraph();
  setGraphAttr(grTE11p);
  grTE11p->SetMarkerStyle(20);
  grTE11p->SetMarkerColor(kBlue);
  grTE11p->SetPoint(grTE11p->GetN(), fTE111 / 1e6, 0.5);
  grTE11p->SetPoint(grTE11p->GetN(), fTE112 / 1e6, 0.5);
  grTE11p->SetPoint(grTE11p->GetN(), fTE113 / 1e6, 0.5);
  grTE11p->SetPoint(grTE11p->GetN(), fTE114 / 1e6, 0.5);
  auto tTE11p = new TText(fTE111 / 1e6, 0.55, "TE11p");
  tTE11p->SetTextColor(kBlue);
  fout.cd();
  tTE11p->Write("tTE11p");

  // TE21p modes
  const double fTE211{cavity.GetResonantModeF(CircularCavity::kTE, 2, 1, 1)};
  const double fTE212{cavity.GetResonantModeF(CircularCavity::kTE, 2, 1, 2)};
  auto grTE21p = new TGraph();
  setGraphAttr(grTE21p);
  grTE21p->SetMarkerStyle(20);
  grTE21p->SetMarkerColor(kBlack);
  grTE21p->SetPoint(grTE21p->GetN(), fTE211 / 1e6, 0.5);
  grTE21p->SetPoint(grTE21p->GetN(), fTE212 / 1e6, 0.5);
  auto tTE21p = new TText(fTE211 / 1e6, 0.55, "TE21p");
  tTE21p->SetTextColor(kBlack);
  fout.cd();
  tTE21p->Write("tTE21p");

  // TM01p modes
  const double fTM010{cavity.GetResonantModeF(CircularCavity::kTM, 0, 1, 0)};
  const double fTM011{cavity.GetResonantModeF(CircularCavity::kTM, 0, 1, 1)};
  const double fTM012{cavity.GetResonantModeF(CircularCavity::kTM, 0, 1, 2)};
  const double fTM013{cavity.GetResonantModeF(CircularCavity::kTM, 0, 1, 3)};
  auto grTM01p = new TGraph();
  setGraphAttr(grTM01p);
  grTM01p->SetMarkerStyle(20);
  grTM01p->SetMarkerColor(kOrange + 1);
  grTM01p->SetPoint(grTM01p->GetN(), fTM010 / 1e6, 0.5);
  grTM01p->SetPoint(grTM01p->GetN(), fTM011 / 1e6, 0.5);
  grTM01p->SetPoint(grTM01p->GetN(), fTM012 / 1e6, 0.5);
  grTM01p->SetPoint(grTM01p->GetN(), fTM013 / 1e6, 0.5);
  auto tTM01p = new TText(fTM010 / 1e6, 0.55, "TM01p");
  tTM01p->SetTextColor(kOrange + 1);
  fout.cd();
  tTM01p->Write("tTM01p");

  // TM11p modes
  const double fTM110{cavity.GetResonantModeF(CircularCavity::kTM, 1, 1, 0)};
  const double fTM111{cavity.GetResonantModeF(CircularCavity::kTM, 1, 1, 1)};
  const double fTM112{cavity.GetResonantModeF(CircularCavity::kTM, 1, 1, 2)};
  const double fTM113{cavity.GetResonantModeF(CircularCavity::kTM, 1, 1, 3)};
  auto grTM11p = new TGraph();
  setGraphAttr(grTM11p);
  grTM11p->SetMarkerStyle(20);
  grTM11p->SetMarkerColor(kMagenta + 1);
  grTM11p->SetPoint(grTM11p->GetN(), fTM110 / 1e6, 1);
  grTM11p->SetPoint(grTM11p->GetN(), fTM111 / 1e6, 1);
  grTM11p->SetPoint(grTM11p->GetN(), fTM112 / 1e6, 1);
  grTM11p->SetPoint(grTM11p->GetN(), fTM113 / 1e6, 1);
  auto tTM11p = new TText(fTM110 / 1e6, 1.05, "TM11p");
  tTM11p->SetTextColor(kMagenta + 1);
  fout.cd();
  tTM11p->Write("tTM11p");

  auto mg = new TMultiGraph("mg", "d = 1.1 m, h = 11 m; f [MHz]; ");
  mg->Add(grTE01p);
  mg->Add(grTE11p);
  mg->Add(grTE21p);
  mg->Add(grTM01p);
  mg->Add(grTM11p);

  fout.cd();
  mg->Write();

  fout.Close();
  return 0;
}