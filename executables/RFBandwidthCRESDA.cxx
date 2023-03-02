/*
  RFBandwidthCRESDA.cxx

  Check the axial frequency of harmonic trap for CRESDA
*/

#include <math.h>

#include <iostream>

#include "BasicFunctions/BasicFunctions.h"
#include "BasicFunctions/Constants.h"
#include "ElectronDynamics/QTNMFields.h"
#include "ElectronDynamics/TrajectoryGen.h"
#include "TF1.h"
#include "TFile.h"
#include "TMath.h"
#include "TString.h"
#include "TTree.h"
#include "TVector3.h"

using namespace rad;
using std::cout;
using std::endl;

int main(int argc, char *argv[]) {
  TString outputFile{argv[1]};
  TFile *fout = new TFile(outputFile, "RECREATE");

  const double bkgB{1};                                    // Tesla
  const double trapB{4e-3};                                // Tesla
  const double coilRadius{0.035};                          // Metres
  const double coilCurrent{2 * trapB * coilRadius / MU0};  // Amps

  HarmonicField *field{new HarmonicField(coilRadius, coilCurrent, bkgB)};
  const TVector3 centralField{
      field->evaluate_field_at_point(TVector3(0, 0, 0))};
  cout << "Central field = " << centralField.Mag() << "T\n";

  TGraph *grB{new TGraph()};
  setGraphAttr(grB);
  grB->GetXaxis()->SetTitle("z [cm]");
  grB->GetYaxis()->SetTitle("|B| [T]");
  for (int n{0}; n < 402; n++) {
    double z{-0.1 + 0.2 * double(n) / double(401)};
    double bMag{field->evaluate_field_magnitude(TVector3(0, 0, z))};
    grB->SetPoint(n, z * 1e2, bMag);
  }
  fout->cd();
  grB->Write("grB");

  // Set up pitch angle scan
  const double trapAngleMin{asin(sqrt(1 - trapB / bkgB))};
  const double pitchAngleStart{1.001 * trapAngleMin};
  const double pitchAngleEnd{90 * TMath::Pi() / 180};
  const int nPnts{2};
  const double simTime{1e-6};
  const double simStepSize{1e-12};
  cout << "Min. pitch angle = " << trapAngleMin * 180 / TMath::Pi()
       << " degrees\n";

  // Electron kinematics
  const double electronKE{18600};  // eV
  const double electronSpeed{GetSpeedFromKE(electronKE, ME)};
  const double tau{2 * R_E / (3 * TMath::C())};

  TGraph *grAxFreq{new TGraph()};
  setGraphAttr(grAxFreq);
  grAxFreq->SetTitle("4 mT harmonic trap; #theta [degrees]; f_{a} [MHz]");
  grAxFreq->SetMarkerStyle(20);

  TGraph *grDeltaF = new TGraph();
  setGraphAttr(grDeltaF);
  grDeltaF->SetMarkerStyle(20);
  grDeltaF->SetTitle("4 mT harmonic trap; #theta [degrees]; #Delta f [Hz]");

  for (int iPnt{0}; iPnt < nPnts; iPnt++) {
    cout << iPnt << endl;
    double thisAngle{pitchAngleStart + (pitchAngleEnd - pitchAngleStart) *
                                           double(iPnt) / double(nPnts - 1)};
    TVector3 V0(electronSpeed * sin(thisAngle), 0,
                electronSpeed * cos(thisAngle));
    const double gyroradius{GetGyroradius(
        V0, field->evaluate_field_at_point(TVector3(0, 0, 0)), ME)};
    TVector3 X0(0, -gyroradius, 0);

    TString trackFile{Form(
        "/home/sjones/work/qtnm/outputs/RFBandwidthCRESDA/track%d.root", iPnt)};
    ElectronTrajectoryGen traj(trackFile, field, X0, V0, simStepSize, simTime,
                               0.0, tau);

    // Now open track file
    TFile *fin = new TFile(trackFile, "READ");
    TTree *tr = (TTree *)fin->Get("tree");
    double time;
    double xPos, yPos, zPos;
    tr->SetBranchAddress("time", &time);
    tr->SetBranchAddress("xPos", &xPos);
    tr->SetBranchAddress("yPos", &yPos);
    tr->SetBranchAddress("zPos", &zPos);

    double bMean{0};
    TGraph *grZ{new TGraph()};
    setGraphAttr(grZ);
    grZ->SetTitle(Form("#theta = %.2f", thisAngle * 180 / TMath::Pi()));
    // Loop over tree entries
    for (int e{0}; e < tr->GetEntries(); e++) {
      tr->GetEntry(e);
      TVector3 ePos(xPos, yPos, zPos);
      bMean += field->evaluate_field_magnitude(ePos);
      grZ->SetPoint(e, time, zPos);
    }
    bMean /= double(tr->GetEntries());

    double f{CalcCyclotronFreq(electronKE, bMean)};
    grDeltaF->SetPoint(iPnt, thisAngle * 180 / TMath::Pi(), f);

    TGraph *grZPgram{MakePowerSpectrumPeriodogram(grZ)};
    double maxFreq{-DBL_MAX};
    double maxPower{-DBL_MAX};
    for (int n{0}; n < grZPgram->GetN(); n++) {
      if (grZPgram->GetPointY(n) > maxPower) {
        maxFreq = grZPgram->GetPointX(n);
        maxPower = grZPgram->GetPointY(n);
      }
    }
    cout << "Axial frequency = " << maxFreq / 1e6 << " MHz\n";
    grAxFreq->SetPoint(iPnt, thisAngle * 180 / TMath::Pi(), maxFreq / 1e6);

    fout->cd();
    grZ->Write(Form("grZ%d", iPnt));
    grZPgram->Write(Form("grZPgram%d", iPnt));

    delete grZ;
    delete grZPgram;

    delete tr;
    fin->Close();
    delete fin;

    cout << "\n";
  }

  fout->cd();
  grAxFreq->Write("grAxFreq");

  for (int n{0}; n < grDeltaF->GetN(); n++) {
    double y{grDeltaF->GetPointY(n)};
    grDeltaF->SetPointY(n, y - grDeltaF->GetPointY(grDeltaF->GetN() - 1));
  }
  grDeltaF->Write("grDeltaF");

  fout->Close();
  delete fout;
  return 0;
}