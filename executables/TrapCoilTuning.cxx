/*
  TrapCoilTuning.cxx

  Look at at some correction coils to make a nice compact trap
*/

#include <iostream>

#include "BasicFunctions/BasicFunctions.h"
#include "BasicFunctions/Constants.h"
#include "ElectronDynamics/QTNMFields.h"
#include "TFile.h"
#include "TGraph.h"
#include "TMultiGraph.h"
#include "TString.h"
#include "TVector3.h"

using namespace rad;

int main() {
  TString outputFile{"~/work/qtnm/MagnetDesigns/TrapCoilTuning/plots.root"};
  TFile fout(outputFile, "recreate");

  const double rCoil{15e-3};         // mm
  const double outerCoilSep{50e-3};  // mm
  const double trapDepth{5.5e-3};    // T
  const double iOuterCoil{2 * trapDepth * rCoil / MU0};
  const double outerCoil1Pos{-outerCoilSep / 2};
  const double outerCoil2Pos{outerCoilSep / 2};
  const double innerCoilOffset{6e-3};
  const double innerCoil1Pos{outerCoil1Pos + innerCoilOffset};
  const double innerCoil2Pos{outerCoil2Pos - innerCoilOffset};
  std::cout << "Coils z pos (mm):\t" << outerCoil1Pos << "\t" << innerCoil1Pos
            << "\t" << innerCoil2Pos << "\t" << outerCoil2Pos << std::endl;

  const double innerTrapDepth{3.1e-3};
  const double iInnerCoil{2 * innerTrapDepth * rCoil / MU0};

  CoilField outerCoil1(rCoil, iOuterCoil, outerCoil1Pos);
  CoilField outerCoil2(rCoil, iOuterCoil, outerCoil2Pos);
  CoilField innerCoil1(rCoil, -iInnerCoil, innerCoil1Pos);
  CoilField innerCoil2(rCoil, -iInnerCoil, innerCoil2Pos);

  FourCoilField fourCoils(rCoil, iOuterCoil, outerCoil2Pos, rCoil, -iInnerCoil,
                          innerCoil2Pos);

  // A central coil
  const double centralTrapDepth{0.3e-3};
  const double iCentralCoil{2 * centralTrapDepth * rCoil / MU0};
  CoilField centralCoil(rCoil, iCentralCoil, 0);

  auto grOuterCoilField = new TGraph();
  setGraphAttr(grOuterCoilField);
  grOuterCoilField->SetTitle("Outer coils; z [mm]; |B| [mT]");
  grOuterCoilField->SetLineWidth(3);
  auto grInnerCoilField = new TGraph();
  setGraphAttr(grInnerCoilField);
  grInnerCoilField->SetTitle("Inner coils; z [mm]; |B| [mT]");
  grInnerCoilField->SetLineWidth(3);
  grInnerCoilField->SetLineColor(kRed);
  auto grCentralCoilField = new TGraph();
  setGraphAttr(grCentralCoilField);
  grCentralCoilField->SetTitle("Central coil; z [mm]; |B| [mT]");
  grCentralCoilField->SetLineWidth(3);
  grCentralCoilField->SetLineColor(kOrange + 1);
  auto gr3CoilField = new TGraph();
  setGraphAttr(gr3CoilField);
  gr3CoilField->SetTitle("All coils; z [mm]; |B| [mT]");
  gr3CoilField->SetLineWidth(3);
  gr3CoilField->SetLineColor(kBlue);
  auto gr4CoilField = new TGraph();
  setGraphAttr(gr4CoilField);
  gr4CoilField->SetTitle("All coils; z [mm]; |B| [mT]");
  gr4CoilField->SetLineWidth(3);
  gr4CoilField->SetLineColor(kBlue);

  const uint nAxPnts{200};
  const double zPosMin{-30e-3};
  const double zPosMax{30e-3};
  for (uint iZ{0}; iZ < nAxPnts; iZ++) {
    const double z{zPosMin +
                   (zPosMax - zPosMin) * double(iZ) / double(nAxPnts - 1)};
    TVector3 pos(0, 0, z);
    TVector3 oCoilField1{outerCoil1.evaluate_field_at_point(pos)};
    TVector3 oCoilField2{outerCoil2.evaluate_field_at_point(pos)};
    TVector3 iCoilField1{innerCoil1.evaluate_field_at_point(pos)};
    TVector3 iCoilField2{innerCoil2.evaluate_field_at_point(pos)};
    TVector3 cCoilField{centralCoil.evaluate_field_at_point(pos)};
    TVector3 fCoilField{fourCoils.evaluate_field_at_point(pos)};
    grOuterCoilField->SetPoint(iZ, z * 1e3,
                               (oCoilField1 + oCoilField2).Mag() * 1e3);
    grInnerCoilField->SetPoint(
        iZ, z * 1e3, (iCoilField1 + iCoilField2 + cCoilField).Mag() * 1e3);
    grCentralCoilField->SetPoint(iZ, z * 1e3, cCoilField.Mag() * 1e3);
    gr3CoilField->SetPoint(
        iZ, z * 1e3, (oCoilField1 + oCoilField2 + cCoilField).Mag() * 1e3);
    gr4CoilField->SetPoint(iZ, z * 1e3, fCoilField.Mag() * 1e3);
  }
  fout.cd();
  grOuterCoilField->Write("grOuterCoilField");
  grInnerCoilField->Write("grInnerCoilField");
  grCentralCoilField->Write("grInnerCoilField");
  gr3CoilField->Write("gr3CoilField");
  gr4CoilField->Write("gr4CoilField");

  auto mg1 = new TMultiGraph("mg1", "3 coil fields; z [mm]; |B| [mT]");
  mg1->Add(grOuterCoilField);
  mg1->Add(grCentralCoilField);
  mg1->Add(gr3CoilField);
  mg1->Write();

  auto mg2 = new TMultiGraph("mg2", "4 coil fields; z [mm]; |B| [mT]");
  mg2->Add(grOuterCoilField);
  mg2->Add(grInnerCoilField);
  mg2->Add(gr4CoilField);
  mg2->Write();

  fout.Close();
  return 0;
}