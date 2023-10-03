/*
  LeeWhitingCoils.cxx

  Explore a 4 trap coil setup

  S. Jones 25/09/23
*/

#include <boost/math/special_functions/legendre.hpp>

#include "BasicFunctions/BasicFunctions.h"
#include "BasicFunctions/Constants.h"
#include "ElectronDynamics/QTNMFields.h"
#include "TCanvas.h"
#include "TF1.h"
#include "TFile.h"
#include "TGraph.h"
#include "TLine.h"
#include "TMultiGraph.h"
#include "TString.h"

using namespace rad;

/// @brief Fit function for a series of Legendre polynomials
/// @param x
/// @param par
/// @return
double BFieldLegendre(double *x, double *par) {
  double xScaled{x[0] / 100};
  double b0{boost::math::legendre_p<double>(0, xScaled) * par[0]};
  double b2{boost::math::legendre_p<double>(2, xScaled) * par[1]};
  double b4{boost::math::legendre_p<double>(4, xScaled) * par[2]};
  double b6{boost::math::legendre_p<double>(6, xScaled) * par[3]};
  double b8{boost::math::legendre_p<double>(8, xScaled) * par[4]};
  double b10{boost::math::legendre_p<double>(4, xScaled) * par[5]};
  return b0 + b2 + b4 + b6 + b8 + b10;
}

int main() {
  TString outputFile{
      "/home/sjones/work/qtnm/MagnetDesigns/LeeWhitingCoils/plots.root"};
  TFile fout(outputFile, "recreate");

  const double lInnerCoil{5e-2};  // m
  const double lOuterCoil{5e-2};  // m

  // Minimum dimensions
  const double minCoilSeparation{204e-3};                          // m
  const double minCoilDiameter{minCoilSeparation / (0.1216 * 2)};  // m

  // Where are our coils positioned?
  const double coilPosition1{0.1216 * minCoilDiameter};
  const double coilPosition2{0.4704 * minCoilDiameter};

  const double centralField{17.96 / minCoilDiameter};  // uT / (A * turns)
  std::cout << centralField << "uT / (A * turns)\n";
  // Therefore for a given field, can calculate the required ampere turns
  const double reqField{0.7};  // T
  const double reqAmpTurns{reqField / (centralField * 1e-6)};
  std::cout << "Required ampere turns = " << reqAmpTurns << std::endl;

  const double n{430};  // turns per metre
  const double turnsInnerCoil{lInnerCoil * n};
  const double turnsOuterCoil{lOuterCoil * n};
  const double iInnerCoil{reqAmpTurns * 4 / turnsInnerCoil};
  const double iOuterCoil{iInnerCoil * 2.2604};

  std::cout << "Coil radii = " << minCoilDiameter * 0.5 * 1e2 << " cm\n";
  std::cout << "Coil edges at:\t" << -coilPosition2 * 1e3 << " mm\t"
            << -coilPosition1 * 1e3 << " mm\t" << coilPosition1 * 1e3 << " mm\t"
            << coilPosition2 * 1e3 << " mm\n";
  std::cout << "Turns per coil:\t" << turnsOuterCoil << "\t" << turnsInnerCoil
            << "\t" << turnsInnerCoil << "\t" << turnsOuterCoil << std::endl;
  std::cout << "Current:\t" << iOuterCoil << " A\t" << iInnerCoil << " A\t"
            << iInnerCoil << " A\t" << iOuterCoil << " A\n";

  // Define the coils
  SolenoidField fieldInner1(minCoilDiameter / 2, lInnerCoil, iInnerCoil, n, MU0,
                            -coilPosition1);
  SolenoidField fieldInner2(minCoilDiameter / 2, lInnerCoil, iInnerCoil, n, MU0,
                            coilPosition1);
  SolenoidField fieldOuter1(minCoilDiameter / 2, lOuterCoil, iOuterCoil, n, MU0,
                            -coilPosition2);
  SolenoidField fieldOuter2(minCoilDiameter / 2, lOuterCoil, iOuterCoil, n, MU0,
                            coilPosition2);

  // Draw the fields
  const uint nPnts{701};
  const uint middlePnt{nPnts / 2};
  std::cout << "Middle point = " << middlePnt << std::endl;
  const double zMin{-0.5};
  const double zMax{0.5};

  auto grInner1 = new TGraph();
  setGraphAttr(grInner1);
  grInner1->SetLineWidth(3);
  grInner1->SetTitle("Inner coil 1; z [mm]; |B| [T]");
  grInner1->SetLineColor(kYellow + 1);
  auto grInner2 = new TGraph();
  setGraphAttr(grInner2);
  grInner2->SetLineWidth(3);
  grInner2->SetTitle("Inner coil 2; z [mm]; |B| [T]");
  grInner2->SetLineColor(kYellow + 1);

  auto grOuter1 = new TGraph();
  setGraphAttr(grOuter1);
  grOuter1->SetLineWidth(3);
  grOuter1->SetTitle("Outer coil 1; z [mm]; |B| [T]");
  grOuter1->SetLineColor(kMagenta + 1);
  auto grOuter2 = new TGraph();
  setGraphAttr(grOuter2);
  grOuter2->SetLineWidth(3);
  grOuter2->SetTitle("Outer coil 2; z [mm]; |B| [T]");
  grOuter2->SetLineColor(kMagenta + 1);

  auto grAllCoils = new TGraph();
  setGraphAttr(grAllCoils);
  grAllCoils->SetLineWidth(3);
  grAllCoils->SetTitle("All coils; z [mm]; |B| [T]");
  grAllCoils->SetLineColor(kRed);

  for (uint iZ{0}; iZ < nPnts; iZ++) {
    const double z{zMin + (zMax - zMin) * double(iZ) / double(nPnts - 1)};
    const double bMagInner1{
        fieldInner1.evaluate_field_magnitude(TVector3(0, 0, z))};
    const double bMagInner2{
        fieldInner2.evaluate_field_magnitude(TVector3(0, 0, z))};
    const double bMagOuter1{
        fieldOuter1.evaluate_field_magnitude(TVector3(0, 0, z))};
    const double bMagOuter2{
        fieldOuter2.evaluate_field_magnitude(TVector3(0, 0, z))};
    const double bMagAll{
        (fieldInner1.evaluate_field_at_point(TVector3(0, 0, z)) +
         fieldInner2.evaluate_field_at_point(TVector3(0, 0, z)) +
         fieldOuter1.evaluate_field_at_point(TVector3(0, 0, z)) +
         fieldOuter2.evaluate_field_at_point(TVector3(0, 0, z)))
            .Mag()};
    grInner1->SetPoint(grInner1->GetN(), z * 1e3, bMagInner1);
    grInner2->SetPoint(grInner2->GetN(), z * 1e3, bMagInner2);
    grOuter1->SetPoint(grOuter1->GetN(), z * 1e3, bMagOuter1);
    grOuter2->SetPoint(grOuter2->GetN(), z * 1e3, bMagOuter2);
    grAllCoils->SetPoint(grAllCoils->GetN(), z * 1e3, bMagAll);

    // if (bMagAll > allCoilsMax) allCoilsMax = bMagAll;
  }
  // Calculate the maximum field
  double allCoilsMax{grAllCoils->GetPointY(middlePnt)};

  fout.cd();
  grInner1->Write("grInner1");
  grInner2->Write("grInner2");
  grOuter1->Write("grOuter1");
  grOuter2->Write("grOuter2");
  grAllCoils->Write("grAllCoils");

  auto mg = new TMultiGraph("mg", Form("R_{coil} = %.0f mm; z [mm]; |B| [T]",
                                       minCoilDiameter * 1e3 / 2));
  mg->Add(grInner1);
  mg->Add(grInner2);
  mg->Add(grOuter1);
  mg->Add(grOuter2);
  mg->Add(grAllCoils);

  auto lInner1_1 = new TLine((-coilPosition1 - lInnerCoil / 2) * 1e3, 0.15,
                             (-coilPosition1 - lInnerCoil / 2) * 1e3, 0.25);
  lInner1_1->SetLineWidth(2);
  lInner1_1->SetLineColor(kYellow + 1);
  lInner1_1->Write("lInner1_1");
  auto lInner1_2 = new TLine((-coilPosition1 + lInnerCoil / 2) * 1e3, 0.15,
                             (-coilPosition1 + lInnerCoil / 2) * 1e3, 0.25);
  lInner1_2->SetLineWidth(2);
  lInner1_2->SetLineColor(kYellow + 1);
  lInner1_2->Write("lInner1_2");
  auto lInner2_1 = new TLine((coilPosition1 - lInnerCoil / 2) * 1e3, 0.15,
                             (coilPosition1 - lInnerCoil / 2) * 1e3, 0.25);
  lInner2_1->SetLineWidth(2);
  lInner2_1->SetLineColor(kYellow + 1);
  lInner2_1->Write("lInner2_1");
  auto lInner2_2 = new TLine((coilPosition1 + lInnerCoil / 2) * 1e3, 0.15,
                             (coilPosition1 + lInnerCoil / 2) * 1e3, 0.25);
  lInner2_2->SetLineWidth(2);
  lInner2_2->SetLineColor(kYellow + 1);
  lInner2_2->Write("lInner2_2");

  auto lOuter1_1 = new TLine((-coilPosition2 - lOuterCoil / 2) * 1e3, 0.4,
                             (-coilPosition2 - lOuterCoil / 2) * 1e3, 0.5);
  lOuter1_1->SetLineWidth(2);
  lOuter1_1->SetLineColor(kMagenta + 1);
  lOuter1_1->Write("lOuter1_1");
  auto lOuter1_2 = new TLine((-coilPosition2 + lOuterCoil / 2) * 1e3, 0.4,
                             (-coilPosition2 + lOuterCoil / 2) * 1e3, 0.5);
  lOuter1_2->SetLineWidth(2);
  lOuter1_2->SetLineColor(kMagenta + 1);
  lOuter1_2->Write("lOuter1_2");
  auto lOuter2_1 = new TLine((coilPosition2 - lOuterCoil / 2) * 1e3, 0.4,
                             (coilPosition2 - lOuterCoil / 2) * 1e3, 0.5);
  lOuter2_1->SetLineWidth(2);
  lOuter2_1->SetLineColor(kMagenta + 1);
  lOuter2_1->Write("lOuter2_1");
  auto lOuter2_2 = new TLine((coilPosition2 + lOuterCoil / 2) * 1e3, 0.4,
                             (coilPosition2 + lOuterCoil / 2) * 1e3, 0.5);
  lOuter2_2->SetLineWidth(2);
  lOuter2_2->SetLineColor(kMagenta + 1);
  lOuter2_2->Write("lOuter2_2");

  fout.cd();
  mg->Write();

  const double fieldDrop{1e-3};
  double homogeneousLen{0};
  for (uint i{middlePnt}; i < grAllCoils->GetN(); i++) {
    if (grAllCoils->GetPointY(i) < allCoilsMax - fieldDrop) {
      homogeneousLen = grAllCoils->GetPointX(i);
      break;
    }
  }
  homogeneousLen *= 2;
  std::cout << "Homogeneous length = " << homogeneousLen << " mm\n";

  auto lUni1 = new TLine(-homogeneousLen / 2, 0.65, -homogeneousLen / 2, 0.75);
  lUni1->SetLineWidth(2);
  lUni1->SetLineColor(kRed);
  lUni1->SetLineStyle(2);
  auto lUni2 = new TLine(homogeneousLen / 2, 0.65, homogeneousLen / 2, 0.75);
  lUni2->SetLineWidth(2);
  lUni2->SetLineColor(kRed);
  lUni2->SetLineStyle(2);
  lUni1->Write("lUni1");
  lUni2->Write("lUni2");

  // Now see if we can make a nice 2D plot of the uniform region
  const double rhoAxMin{0};
  const double rhoAxMax{0.3};
  const double zAxMin{-0.3};
  const double zAxMax{0.3};
  uint nAxPnts{300};
  auto h2Uniform =
      new TH2D("h2Uniform",
               Form("R_{coil} = %.0f mm; z [mm]; #rho [mm]; |B| [T]",
                    minCoilDiameter * 1e3 / 2),
               nAxPnts, zAxMin * 1e3, zAxMax * 1e3, nAxPnts, rhoAxMin * 1e3,
               rhoAxMax * 1e3);
  SetHistAttr(h2Uniform);
  auto h2FieldDiff =
      new TH2D("h2FieldDiff",
               Form("R_{coil} = %.0f mm; z [mm]; #rho [mm]; #DeltaB / B",
                    minCoilDiameter * 1e3 / 2),
               nAxPnts, zAxMin * 1e3, zAxMax * 1e3, nAxPnts, rhoAxMin * 1e3,
               rhoAxMax * 1e3);
  SetHistAttr(h2FieldDiff);

  for (uint iZ{1}; iZ <= h2Uniform->GetNbinsX(); iZ++) {
    double z{h2Uniform->GetXaxis()->GetBinCenter(iZ)};  // mm
    for (uint iY{1}; iY <= h2Uniform->GetNbinsY(); iY++) {
      double y{h2Uniform->GetYaxis()->GetBinCenter(iY)};  // mm
      TVector3 pos(0, y / 1e3, z / 1e3);
      TVector3 bField{fieldInner1.evaluate_field_at_point(pos) +
                      fieldInner2.evaluate_field_at_point(pos) +
                      fieldOuter1.evaluate_field_at_point(pos) +
                      fieldOuter2.evaluate_field_at_point(pos)};

      h2Uniform->SetBinContent(iZ, iY, bField.Mag());
      h2FieldDiff->SetBinContent(iZ, iY,
                                 abs(allCoilsMax - bField.Mag()) / allCoilsMax);
    }
  }

  double contours[5];
  contours[0] = allCoilsMax - 1e-3;
  contours[1] = allCoilsMax * (1.0 - 100e-6);
  contours[2] = allCoilsMax * (1.0 - 10e-6);
  contours[3] = allCoilsMax * (1.0 - 1e-6);
  contours[4] = allCoilsMax;
  // contours[1] = allCoilsMax + 1e-3;
  h2Uniform->SetContour(5, contours);
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

  /*
  const double fitLimits{100};  // mm
  auto fLeg = new TF1("fLeg", BFieldLegendre, -fitLimits, fitLimits, 6);
  fLeg->SetParameter(0, 0.7);
  fLeg->SetParameter(1, -1);
  fLeg->SetParameter(2, -0.001);
  fLeg->SetParameter(3, -0.001);
  fLeg->SetParameter(4, -0.001);
  fLeg->SetParameter(5, -0.001);
  grAllCoils->Fit(fLeg, "r");

  fLeg->Write();
  */

  std::cout << "Central field max = " << allCoilsMax << " T\n";

  fout.Close();
  return 0;
}