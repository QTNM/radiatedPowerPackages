/*
  HelmholtzCoils.cxx

  See what kind of homogeneous region we can can from a system of Helholtz coils
*/

#include <iostream>

#include "BasicFunctions/BasicFunctions.h"
#include "BasicFunctions/Constants.h"
#include "ElectronDynamics/QTNMFields.h"
#include "TFile.h"
#include "TGraph.h"
#include "TH2.h"
#include "TLine.h"
#include "TMultiGraph.h"
#include "TString.h"

using namespace rad;

int main() {
  TString outputFile{
      "/home/sjones/work/qtnm/MagnetDesigns/HelmholtzCoils/fieldPlots.root"};
  TFile fout(outputFile, "recreate");

  const double bField{0.7 / (2 * pow(0.8, 1.5))};  // Tesla
  const double n{430};                             // turns per metre

  const double flangeDiameter{0.11};              // m
  const double lSolenoid{5e-2};                   // m
  const double coilSeparation{204e-3};            // m
  const double rSolenoid{coilSeparation / 1.02};  // m

  const double offset{coilSeparation / 2};  // m
  const double iSolenoid{2 * bField *
                         sqrt(pow(lSolenoid / 2, 2) + rSolenoid * rSolenoid) /
                         (MU0 * n * lSolenoid)};  // Amps
  std::cout << "Coil radii, current = " << rSolenoid << " mm,\t" << iSolenoid
            << " A\n";

  SolenoidField field1(rSolenoid, lSolenoid, iSolenoid, n, MU0, offset);
  SolenoidField field2(rSolenoid, lSolenoid, iSolenoid, n, MU0, -offset);

  const double turnsPerCoil{n * lSolenoid};
  const double centralHelmholtzField{MU0 * turnsPerCoil * iSolenoid /
                                     (2 * rSolenoid) * 1.43108};
  std::cout << "Central Helmholtz field = " << centralHelmholtzField << " T\n";

  const uint nPnts{401};
  const double zMin{-0.2};
  const double zMax{0.2};

  auto grCoil1 = new TGraph();
  setGraphAttr(grCoil1);
  grCoil1->SetLineWidth(3);
  grCoil1->SetTitle("Coil 1; z [mm]; |B| [T]");
  grCoil1->SetLineColor(kYellow + 1);
  auto grCoil2 = new TGraph();
  setGraphAttr(grCoil2);
  grCoil2->SetLineWidth(3);
  grCoil2->SetTitle("Coil 2; z [mm]; |B| [T]");
  grCoil2->SetLineColor(kMagenta + 1);
  auto grBothCoils = new TGraph();
  setGraphAttr(grBothCoils);
  grBothCoils->SetLineWidth(3);
  grBothCoils->SetTitle("Both coils; z [mm]; |B| [T]");
  grBothCoils->SetLineColor(kRed);

  for (uint iZ{0}; iZ < nPnts; iZ++) {
    const double z{zMin + (zMax - zMin) * double(iZ) / double(nPnts - 1)};
    const double bMag1{field1.evaluate_field_magnitude(TVector3(0, 0, z))};
    const double bMag2{field2.evaluate_field_magnitude(TVector3(0, 0, z))};
    const double bMagBoth{(field1.evaluate_field_at_point(TVector3(0, 0, z)) +
                           field2.evaluate_field_at_point(TVector3(0, 0, z)))
                              .Mag()};
    grCoil1->SetPoint(grCoil1->GetN(), z * 1e3, bMag1);
    grCoil2->SetPoint(grCoil2->GetN(), z * 1e3, bMag2);
    grBothCoils->SetPoint(grBothCoils->GetN(), z * 1e3, bMagBoth);
  }
  fout.cd();
  grCoil1->Write("grCoil1");
  grCoil2->Write("grCoil2");
  grBothCoils->Write("grBothCoils");

  auto mg = new TMultiGraph(
      "mg", Form("R_{coil} = %.0f mm; z [mm]; |B| [T]", rSolenoid * 1e3));
  mg->Add(grCoil1);
  mg->Add(grCoil2);
  mg->Add(grBothCoils);
  fout.cd();
  mg->Write();

  // Now get uniform region
  const double maxField{grBothCoils->GetPointY(200)};
  std::cout << "Max field = " << maxField << " T\n";

  double homogeneousLength{0};
  for (uint i{200}; i < grBothCoils->GetN(); i++) {
    if (grBothCoils->GetPointY(i) < maxField - 1e-3) {
      homogeneousLength = grBothCoils->GetPointX(i);
      break;
    }
  }
  homogeneousLength *= 2;
  std::cout << "Homogeneous length = " << homogeneousLength << " mm\n";

  auto l1_1 = new TLine((offset + lSolenoid / 2) * 1e3, 0.8 * bField,
                        (offset + lSolenoid / 2) * 1e3, bField);
  l1_1->SetLineWidth(3);
  l1_1->SetLineColor(kYellow + 1);
  auto l1_2 = new TLine((offset - lSolenoid / 2) * 1e3, 0.8 * bField,
                        (offset - lSolenoid / 2) * 1e3, bField);
  l1_2->SetLineWidth(3);
  l1_2->SetLineColor(kYellow + 1);
  auto l2_1 = new TLine((-offset + lSolenoid / 2) * 1e3, 0.8 * bField,
                        (-offset + lSolenoid / 2) * 1e3, bField);
  l2_1->SetLineWidth(3);
  l2_1->SetLineColor(kMagenta + 1);
  auto l2_2 = new TLine((-offset - lSolenoid / 2) * 1e3, 0.8 * bField,
                        (-offset - lSolenoid / 2) * 1e3, bField);
  l2_2->SetLineWidth(3);
  l2_2->SetLineColor(kMagenta + 1);

  auto lUni1 =
      new TLine(-homogeneousLength / 2, 0.6, -homogeneousLength / 2, 0.8);
  lUni1->SetLineWidth(3);
  lUni1->SetLineStyle(2);
  lUni1->SetLineColor(kRed);
  auto lUni2 =
      new TLine(homogeneousLength / 2, 0.6, homogeneousLength / 2, 0.8);
  lUni2->SetLineWidth(3);
  lUni2->SetLineStyle(2);
  lUni2->SetLineColor(kRed);

  fout.cd();
  l1_1->Write("l1_1");
  l1_2->Write("l1_2");
  l2_1->Write("l2_1");
  l2_2->Write("l2_2");
  lUni1->Write("lUni1");
  lUni2->Write("lUni2");

  // Now try and make a countour plot showing where the uniform region is
  const double rhoAxMin{0};
  const double rhoAxMax{0.1};
  const double zAxMin{-0.1};
  const double zAxMax{0.1};
  uint nAxPnts{300};
  auto h2Uniform = new TH2D(
      "h2Uniform",
      Form("R_{coil} = %.0f mm; z [mm]; #rho [mm]; |B| [T]", rSolenoid * 1e3),
      nAxPnts, zAxMin * 1e3, zAxMax * 1e3, nAxPnts, rhoAxMin * 1e3,
      rhoAxMax * 1e3);
  SetHistAttr(h2Uniform);
  auto h2FieldDiff =
      new TH2D("h2FieldDiff",
               Form("R_{coil} = %.0f mm; z [mm]; #rho [mm]; #DeltaB / B",
                    rSolenoid * 1e3),
               nAxPnts, zAxMin * 1e3, zAxMax * 1e3, nAxPnts, rhoAxMin * 1e3,
               rhoAxMax * 1e3);
  SetHistAttr(h2FieldDiff);

  for (uint iZ{1}; iZ <= h2Uniform->GetNbinsX(); iZ++) {
    double z{h2Uniform->GetXaxis()->GetBinCenter(iZ)};  // mm
    for (uint iY{1}; iY <= h2Uniform->GetNbinsY(); iY++) {
      double y{h2Uniform->GetYaxis()->GetBinCenter(iY)};  // mm
      TVector3 pos(0, y / 1e3, z / 1e3);
      TVector3 bField{field1.evaluate_field_at_point(pos) +
                      field2.evaluate_field_at_point(pos)};

      h2Uniform->SetBinContent(iZ, iY, bField.Mag());
      h2FieldDiff->SetBinContent(iZ, iY,
                                 abs(maxField - bField.Mag()) / maxField);
    }
  }

  double contours[2];
  contours[0] = maxField - 1e-3;
  // contours[2] = maxField;
  contours[1] = maxField + 1e-3;
  h2Uniform->SetContour(2, contours);
  h2Uniform->SetLineColor(kBlack);
  h2Uniform->SetLineWidth(3);
  h2Uniform->Write();

  double contoursDiff[5];
  contoursDiff[0] = 1e-6;
  contoursDiff[1] = 10e-6;
  contoursDiff[2] = 100e-6;
  contoursDiff[3] = 1000e-6;
  contoursDiff[4] = 0.0014;
  h2FieldDiff->SetContour(5, contoursDiff);
  h2FieldDiff->Write();

  fout.Close();
  return 0;
}