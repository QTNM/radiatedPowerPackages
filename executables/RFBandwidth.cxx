/*
  RFBandwdth.cxx

  Check the required RF bandwidth for a hypothetical metre scale experiment

  Command line argument is just output ROOT file
*/

#include <math.h>

#include <iostream>

#include "BasicFunctions/BasicFunctions.h"
#include "BasicFunctions/Constants.h"
#include "ElectronDynamics/QTNMFields.h"
#include "ElectronDynamics/TrajectoryGen.h"
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

  const double trapFraction{0.1};  //  Desired fraction of electrons to trap
  const double dThetaMax{asin(trapFraction)};
  const double trapAngleMin{TMath::Pi() / 2 - dThetaMax};
  cout << "dThetaMax = " << dThetaMax * 180 / TMath::Pi() << "degrees\n";
  cout << "Min. pitch angle = " << trapAngleMin * 180 / TMath::Pi()
       << " degrees\n";

  const double bkgBMin{1};
  const double bkgBMax{1.0};
  const int nbkgBPnts{2};

  for (int iBkg{0}; iBkg < nbkgBPnts; iBkg++) {
    // Background field in T
    double bkgB{bkgBMin +
                (bkgBMax - bkgBMin) * double(iBkg) / double(nbkgBPnts - 1)};
    const double trapB{bkgB / pow(cos(dThetaMax), 2) -
                       bkgB};  // Trapping field in T
    cout << "Background field = " << bkgB
         << " T,\t Trap depth = " << trapB * 1e3 << " mT\n";

    const double coilRadius{0.1};                            // Metres
    const double coilCurrent{2 * trapB * coilRadius / MU0};  // Amps
    const double trapLength{1.0};                            // Metres

    BathtubField *field1m =
        new BathtubField(coilRadius, coilCurrent, -trapLength / 2,
                         trapLength / 2, TVector3(0, 0, bkgB));
    const TVector3 centralField{
        field1m->evaluate_field_at_point(TVector3(0, 0, 0))};
    cout << "Central field = " << centralField.Mag() << "T\n";

    // Draw field
    TGraph *grB{new TGraph()};
    setGraphAttr(grB);
    grB->GetXaxis()->SetTitle("z [m]");
    grB->GetYaxis()->SetTitle("|B| [T]");
    for (int n{0}; n < 402; n++) {
      double z{-0.55 + 1.1 * double(n) / double(401)};
      double bMag{field1m->evaluate_field_magnitude(TVector3(0, 0, z))};
      grB->SetPoint(n, z, bMag);
    }
    fout->cd();
    grB->Write("grB");

    const double pitchAngleStart{1.001 * trapAngleMin};
    const double pitchAngleEnd{90 * TMath::Pi() / 180};
    const int nPnts{20};
    const double simTime{5e-6};
    const double simStepSize{1e-12};

    // Electron kinematics
    const double electronKE{18600};  // eV
    const double electronSpeed{GetSpeedFromKE(electronKE, ME)};
    const double tau = 2 * R_E / (3 * TMath::C());

    TGraph *grBMean = new TGraph();
    setGraphAttr(grBMean);
    grBMean->SetMarkerStyle(20);
    grBMean->SetTitle(
        Form("1 m long bathtub trap, B_{bkg} = %.1f T, #Delta B = %.1f mT; "
             "#theta [degrees]; B_{mean} [T]",
             bkgB, trapB * 1e3));

    TGraph *grF = new TGraph();
    setGraphAttr(grF);
    grF->SetMarkerStyle(20);
    grF->SetTitle(
        Form("1 m long bathtub trap, B_{bkg} = %.1f T, #Delta B = %.1f mT; "
             "#theta [degrees]; f [Hz]",
             bkgB, trapB * 1e3));

    TGraph *grDeltaF = new TGraph();
    setGraphAttr(grDeltaF);
    grDeltaF->SetMarkerStyle(20);
    grDeltaF->SetTitle(Form(
        "B_{bkg} = %.1f T, #Delta B = %.1f mT; #theta [degrees]; #Delta f [Hz]",
        bkgB, trapB * 1e3));
    grDeltaF->SetMarkerColor(52 + iBkg * 3);

    TGraph *grAxFreq = new TGraph();
    setGraphAttr(grAxFreq);
    grAxFreq->SetMarkerStyle(20);
    grAxFreq->SetTitle(Form(
        "B_{bkg} = %.1f T, #Delta B = %.1f mT; #theta [degrees]; f_{a} [MHz]",
        bkgB, trapB * 1e3));

    for (int iPnt{0}; iPnt < nPnts; iPnt++) {
      double thisAngle{pitchAngleStart + (pitchAngleEnd - pitchAngleStart) *
                                             double(iPnt) / double(nPnts - 1)};
      TVector3 V0(electronSpeed * sin(thisAngle), 0,
                  electronSpeed * cos(thisAngle));
      const double gyroradius{GetGyroradius(
          V0, field1m->evaluate_field_at_point(TVector3(0, 0, 0)), ME)};
      TVector3 X0(0, -gyroradius, 0);
      TString trackFile{Form(
          "/home/sjones/work/qtnm/outputs/RFBandwidth/track%d.root", iPnt)};
      ElectronTrajectoryGen traj(trackFile, field1m, X0, V0, simStepSize,
                                 simTime);

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
        bMean += field1m->evaluate_field_magnitude(ePos);
        grZ->SetPoint(e, time, zPos);
      }
      bMean /= double(tr->GetEntries());
      cout << iPnt << ":\t bMean = " << bMean << " T\n";

      double f{CalcCyclotronFreq(electronKE, bMean)};
      grBMean->SetPoint(iPnt, thisAngle * 180 / TMath::Pi(), bMean);
      grF->SetPoint(iPnt, thisAngle * 180 / TMath::Pi(), f);
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
      cout << "Axial frequency = " << maxFreq / 1e6 << "MHz\n";
      grAxFreq->SetPoint(iPnt, thisAngle * 180 / TMath::Pi(), maxFreq / 1e6);

      fout->cd();
      grZ->Write(Form("grZ_%d_%d", iBkg, iPnt));
      delete grZ;
      delete grZPgram;

      delete tr;
      fin->Close();
      delete fin;
    }

    for (int n{0}; n < grDeltaF->GetN(); n++) {
      double y{grDeltaF->GetPointY(n)};
      grDeltaF->SetPointY(n, y - grDeltaF->GetPointY(grDeltaF->GetN() - 1));
    }

    fout->cd();
    grBMean->Write(Form("grBMean%d", iBkg));
    grF->Write(Form("grF%d", iBkg));
    grDeltaF->Write(Form("grDeltaF%d", iBkg));
    grAxFreq->Write(Form("grAxFreq%d", iBkg));

    delete grBMean;
    delete grF;
    delete grDeltaF;
    cout << "\n";
  }  // Loop over background fields

  fout->Close();
  delete fout;
  return 0;
}