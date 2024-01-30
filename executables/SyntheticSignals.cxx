/*
  SyntheticSignals.cxx

  Comparison of magnetic dipole and crossed electric dipoles to our actual
  signal
*/

#include <getopt.h>

#include <cmath>
#include <iostream>
#include <string>

#include "Antennas/IsotropicAntenna.h"
#include "BasicFunctions/BasicFunctions.h"
#include "BasicFunctions/Constants.h"
#include "ElectronDynamics/QTNMFields.h"
#include "ElectronDynamics/TrajectoryGen.h"
#include "FieldClasses/FieldClasses.h"
#include "TCanvas.h"
#include "TFile.h"
#include "TH2.h"
#include "TString.h"
#include "TStyle.h"
#include "TTree.h"
#include "TVector3.h"

using namespace rad;

TVector3 EFieldMagneticDipole(double r, double theta, double phi, double t,
                              double f, double phase = 0, double m0 = 1) {
  const double omega{2 * M_PI * f};
  const double prefac{MU0 * m0 * sin(theta) / (4 * M_PI * r)};
  const double ePhi{(omega * omega / TMath::C()) *
                        cos(omega * (t - r / TMath::C())) +
                    (omega / r) * sin(omega * (t - r / TMath::C()))};
  return prefac * TVector3(-ePhi * sin(phi), ePhi * cos(phi), 0);
}

TVector3 EFieldMagneticDipole(TVector3 vec, double t, double f,
                              double phase = 0, double m0 = 1) {
  return EFieldMagneticDipole(vec.Mag(), vec.Theta(), vec.Phi(), t, f, phase,
                              m0);
}

TVector3 EFieldElectricDipole(TVector3 r, TVector3 p0, double f, double t,
                              double phase = 0) {
  double omega{2 * M_PI * f};
  TVector3 eUnit{(p0.Cross(r.Unit())).Cross(r.Unit())};
  double prefac{-MU0 * omega * omega / (4 * M_PI * r.Mag())};
  return eUnit * prefac * cos(omega * (t - r.Mag() / TMath::C()) + phase);
}

TVector3 EFieldCrossedDipole(TVector3 r, double f, double t) {
  TVector3 dipole1(1, 0, 0);
  TVector3 dipole2(0, 1, 0);
  double omega{2 * M_PI * f};
  return EFieldElectricDipole(r, dipole1, f, t, 0) +
         EFieldElectricDipole(r, dipole2, f, t, M_PI_2);
}

double PRadMagneticDipole(double f, double m0) {
  double omega{2 * M_PI * f};
  return MU0 * m0 * m0 * pow(omega, 4) / (12 * M_PI * pow(TMath::C(), 3));
}

int main(int argc, char *argv[]) {
  int opt{};
  std::string outputFileStr{" "};
  while ((opt = getopt(argc, argv, ":o:")) != -1) {
    switch (opt) {
      case 'o':
        outputFileStr = optarg;
        break;

      case ':':
        std::cout << "Option needs a value\n";
        break;

      case '?':
        std::cout << "Unknown option: " << optopt << std::endl;
        break;
    }
  }

  TString outputFile{outputFileStr};
  TFile fout(outputFile, "recreate");

  const double bField{0.65};                     // Tesla
  const double eKE{18.6e3};                      // eV
  const double eSpeed{GetSpeedFromKE(eKE, ME)};  // m s^-1
  TVector3 vel(eSpeed, 0, 0);
  const double fCyc{CalcCyclotronFreq(eKE, bField)};        // Hertz
  const double PRad{CalcLarmorPower(eKE, bField, M_PI_2)};  // Watts
  std::cout << "f_cyc = " << fCyc / 1e9 << " GHz\tP_rad = " << PRad * 1e15
            << " fW\n";

  auto field = new UniformField(bField);
  const double r_g{GetGyroradius(vel, TVector3(0, 0, bField), ME)};
  TVector3 x0(0, -r_g, 0);
  const double simTime{1e-7};       // seconds
  const double simStepSize{2e-12};  // seconds
  TString trackFile{"~/work/qtnm/SyntheticSignals/track.root"};
  ElectronTrajectoryGen traj(trackFile, field, x0, vel, simStepSize, simTime);

  // Define source point
  TVector3 fieldPoint(0, 10e-2, 0);
  auto ant = new IsotropicAntenna(fieldPoint, 1, 1, fCyc);
  FieldPoint fp(trackFile, ant);
  fp.GenerateFields(0, simTime);

  auto grRealX{fp.GetEFieldTimeDomain(FieldPoint::Coord_t::kX, true)};
  auto grRealY{fp.GetEFieldTimeDomain(FieldPoint::Coord_t::kY, true)};
  auto grRealZ{fp.GetEFieldTimeDomain(FieldPoint::Coord_t::kZ, true)};

  fout.cd();
  grRealX->Write("grRealX");
  grRealY->Write("grRealY");
  grRealZ->Write("grRealZ");

  // Now try magnetic dipole
  const double kCyc{2 * M_PI * fCyc / TMath::C()};
  const double i0{sqrt((12 * pow(TMath::C(), 3) * PRad) /
                       (MU0 * M_PI * pow(r_g, 4) * pow(2 * M_PI * fCyc, 4)))};
  const double m0{M_PI * r_g * r_g * i0};
  const double magDipoleP{PRadMagneticDipole(fCyc, m0)};
  std::cout << "Magnetic dipole radiated power = " << magDipoleP * 1e15
            << " fW\n";

  // Load in the file
  TFile fin(trackFile, "read");
  auto tr = (TTree *)fin.Get("tree");
  double time{};
  tr->SetBranchAddress("time", &time);
  auto grMagDX = new TGraph();
  setGraphAttr(grMagDX);
  auto grMagDY = new TGraph();
  setGraphAttr(grMagDY);
  auto grMagDZ = new TGraph();
  setGraphAttr(grMagDZ);

  auto grEleCDX = new TGraph();
  setGraphAttr(grEleCDX);
  auto grEleCDY = new TGraph();
  setGraphAttr(grEleCDY);
  auto grEleCDZ = new TGraph();
  setGraphAttr(grEleCDZ);

  for (int iE{0}; iE < tr->GetEntries(); iE++) {
    tr->GetEntry(iE);
    TVector3 E{EFieldMagneticDipole(fieldPoint, time, fCyc, 0, m0)};
    TVector3 ECD{EFieldCrossedDipole(fieldPoint, fCyc, time)};
    grMagDX->SetPoint(grMagDX->GetN(), time, E.X());
    grMagDY->SetPoint(grMagDY->GetN(), time, E.Y());
    grMagDZ->SetPoint(grMagDZ->GetN(), time, E.Z());

    grEleCDX->SetPoint(grEleCDX->GetN(), time, ECD.X());
    grEleCDY->SetPoint(grEleCDY->GetN(), time, ECD.Y());
    grEleCDZ->SetPoint(grEleCDZ->GetN(), time, ECD.Z());
  }
  fin.Close();

  fout.cd();
  grMagDX->Write("grMagDX");
  grMagDY->Write("grMagDY");
  grMagDZ->Write("grMagDZ");

  grEleCDX->Write("grEleCDX");
  grEleCDY->Write("grEleCDY");
  grEleCDZ->Write("grEleCDZ");

  const uint nTimeSteps{40};
  const double timePeriod{1 / fCyc};

  double maxMagD{-DBL_MAX};
  double maxEleCD{-DBL_MAX};
  std::vector<TH2D> MagDVec;
  std::vector<TH2D> EleCDVec;
  for (uint iStep{0}; iStep < nTimeSteps; iStep++) {
    const double thisTime{timePeriod * double(iStep) / double(nTimeSteps)};

    // Now make some 2D plots
    const double nBinsX{200};
    const double nBinsY{150};
    const double xMax{70};
    const double yMax{70};
    TH2D h2MagD(Form("h2MagD_%d", iStep),
                Form("Magnetic dipole: %.1f ps; x [mm]; y [mm]; E_{#phi}",
                     thisTime * 1e12),
                nBinsX, -xMax, xMax, nBinsY, -yMax, yMax);
    SetHistAttr(h2MagD);
    TH2D h2EleCD(
        Form("h2EleCD_%d", iStep),
        Form("Crossed electric dipoles: %.1f ps; x [mm]; y [mm]; E_{#phi}",
             thisTime * 1e12),
        nBinsX, -xMax, xMax, nBinsY, -yMax, yMax);
    SetHistAttr(h2EleCD);
    TH2D h2Real(Form("h2Real_%d", iStep),
                Form("Real electron: %.1f ps; x [mm]; y [mm]; E_{#phi}",
                     thisTime * 1e12),
                nBinsX, -xMax, xMax, nBinsY, -yMax, yMax);
    SetHistAttr(h2Real);

    for (int ix{1}; ix <= h2MagD.GetNbinsX(); ix++) {
      double x{h2MagD.GetXaxis()->GetBinCenter(ix) * 1e-3};
      for (int iy{1}; iy <= h2MagD.GetNbinsY(); iy++) {
        double y{h2MagD.GetYaxis()->GetBinCenter(iy) * 1e-3};
        TVector3 pos(x, y, 0);
        if (pos.Perp() < 5e-3) continue;

        double phi{atan2(y, x)};
        TVector3 phiHat(-sin(phi), cos(phi), 0);

        TVector3 eFieldMagD{EFieldMagneticDipole(pos, thisTime, fCyc, 0, m0)};
        double eFieldMagDDot{eFieldMagD.Dot(phiHat)};
        h2MagD.SetBinContent(ix, iy, eFieldMagDDot);
        if (eFieldMagDDot > maxMagD) maxMagD = eFieldMagDDot;
        TVector3 eFieldEleCD{EFieldCrossedDipole(pos, fCyc, thisTime)};
        double eFieldEleCDDot{eFieldEleCD.Dot(phiHat)};
        if (eFieldEleCDDot > maxEleCD) maxEleCD = eFieldEleCDDot;
        h2EleCD.SetBinContent(ix, iy, eFieldEleCDDot);
      }
    }

    MagDVec.push_back(h2MagD);
    EleCDVec.push_back(h2EleCD);
  }

  gStyle->SetOptStat(0);
  gStyle->SetPalette(56);
  gStyle->SetNumberContours(80);
  for (int ii{0}; ii < MagDVec.size(); ii++) {
    MagDVec.at(ii).Scale(1 / maxMagD);
    MagDVec.at(ii).GetZaxis()->SetRangeUser(-1, 1);
    MagDVec.at(ii).GetZaxis()->SetTitleSize(0.05);
    TCanvas c1(Form("cMagD%d", ii), Form("cMagD%d", ii), 500, 500);
    c1.SetRightMargin(0.18);
    c1.SetLeftMargin(0.15);
    MagDVec.at(ii).Draw("colz");
    c1.Print(Form("~/work/qtnm/SyntheticSignals/cMagD%02d.pdf", ii));
    c1.Print(Form("~/work/qtnm/SyntheticSignals/cMagD%02d.png", ii));
    c1.Print(Form("~/work/qtnm/SyntheticSignals/cMagD%02d.tex", ii));

    EleCDVec.at(ii).Scale(1 / maxEleCD);
    EleCDVec.at(ii).GetZaxis()->SetRangeUser(-1, 1);
    EleCDVec.at(ii).GetZaxis()->SetTitleSize(0.05);
    TCanvas c2(Form("cEleCD%d", ii), Form("cEleCD%d", ii), 500, 500);
    c2.SetRightMargin(0.18);
    c2.SetLeftMargin(0.15);
    EleCDVec.at(ii).Draw("colz");
    c2.Print(Form("~/work/qtnm/SyntheticSignals/cEleCD%02d.pdf", ii));
    c2.Print(Form("~/work/qtnm/SyntheticSignals/cEleCD%02d.png", ii));
    c2.Print(Form("~/work/qtnm/SyntheticSignals/cEleCD%02d.tex", ii));
  }

  fout.Close();
  return 0;
}