/*
  CavityFieldAmps.cxx

  Plot the electric field amplitudes for a given cavity
*/

#include <iostream>

#include "BasicFunctions/BasicFunctions.h"
#include "BasicFunctions/Constants.h"
#include "ElectronDynamics/QTNMFields.h"
#include "ElectronDynamics/TrajectoryGen.h"
#include "TFile.h"
#include "TGraph.h"
#include "TMath.h"
#include "TSpline.h"
#include "TString.h"
#include "TTreeReader.h"
#include "TTreeReaderValue.h"
#include "Waveguides/CircularCavity.h"

using namespace rad;

int main() {
  TString outStr{"/home/sjones/work/qtnm/Cavities/CavityFieldAmps/plots.root"};
  TFile outPlots(outStr, "recreate");

  // Create a harmonic trap
  const double bkgFieldMag{0.7};  // T
  const double rCoil{0.01};       // metres
  const double trapDepth{4e-3};   // T
  const double iCoil{2 * trapDepth * rCoil / MU0};
  auto field = new HarmonicField(rCoil, iCoil, bkgFieldMag);
  // Get the field at the centre of the trap
  TVector3 centreField{field->evaluate_field_at_point(TVector3(0, 0, 0))};

  // Electron kinematics
  const double eKE{18.6e3};  // eV
  const double eSpeed{GetSpeedFromKE(eKE, ME)};

  // Calculate cyclotron frequency
  const double cycFreq{CalcCyclotronFreq(eKE, centreField.Mag())};
  const double k0{TMath::TwoPi() * cycFreq / TMath::C()};
  std::cout << "Cyclotron frequency = " << cycFreq / 1e9 << " GHz\n";

  // Define some cavity stuff
  const double cavityQ{200};
  // Due to non-infinite Q, resonant frequency shifted slightly
  // Calculate what the unshifted resonant frequency is
  const double modeFreq{2 * cavityQ * cycFreq / (2 * cavityQ - 1)};
  std::cout << "Mode resonant frequency = " << modeFreq / 1e9 << " GHz\n";
  const double cavityRadius{5e-3};  // metres
  const double p11Prime{GetBesselPrimeZero(1, 1)};
  const double cavityLength{
      TMath::Pi() / sqrt(pow(TMath::TwoPi() * modeFreq / TMath::C(), 2) -
                         pow(p11Prime / cavityRadius, 2))};
  std::cout << "Cavity length = " << cavityLength * 1e3 << " mm\n";
  // Define the actual cavity
  CircularCavity cav(cavityRadius, cavityLength);
  const double fTE111{cav.GetResonantModeF(CircularCavity::kTE, 1, 1, 1)};
  const double kTE111{TMath::TwoPi() * fTE111 / TMath::C()};
  const double tTE111{1 / fTE111};
  std::cout << "TE111 frequency = " << fTE111 / 1e9 << " GHz\n";

  // Try and word out a normalisation for the TE111 mode
  const uint nNormPnts{30};
  const double dRho{cav.GetRadius() / double(nNormPnts)};
  const double dPhi{TMath::TwoPi() / double(nNormPnts)};
  const double dZ{cav.GetLength() / double(nNormPnts)};
  double integral1{0};
  double integral2{0};
  const double A{1};
  for (uint iRho{0}; iRho < nNormPnts; iRho++) {
    const double rho{dRho / 2 + double(iRho) * dRho};
    for (uint iPhi{0}; iPhi < nNormPnts; iPhi++) {
      const double phi{dPhi / 2 + double(iPhi) * dPhi};
      for (uint iZ{0}; iZ < nNormPnts; iZ++) {
        const double z{-cav.GetLength() / 2 + dZ / 2 + double(iZ) * dZ};
        const double dV{rho * dRho * dPhi * dZ};
        ComplexVector3 eField1{cav.GetModeEField(
            rho, phi, z, CircularCavity::kTE, A, 1, 1, 1, true, 0)};
        integral1 += (eField1.Dot(eField1.Conj())).real() * dV;
        ComplexVector3 eField2{cav.GetModeEField(
            rho, phi, z, CircularCavity::kTE, A, 1, 1, 1, false, tTE111 / 4)};
        integral2 += (eField2.Dot(eField2.Conj())).real() * dV;
      }
    }
  }
  std::cout << "Integral 1, 2 = " << integral1 << ", " << integral2
            << std::endl;
  const double normalisation1{1 / sqrt(integral1)};
  const double normalisation2{1 / sqrt(integral2)};

  // Now we have the normalisations we can propagate our electron
  // Start with an electron completing no axial motion
  const double pitchAngleDeg{89};
  const double pitchAngleRad{pitchAngleDeg * TMath::Pi() / 180};
  TVector3 eVel(eSpeed * sin(pitchAngleRad), 0, eSpeed * cos(pitchAngleRad));
  const double gyroradius{GetGyroradius(
      eVel, field->evaluate_field_at_point(TVector3(0, 0, 0)), ME)};
  std::cout << "Gyroradius = " << gyroradius * 1e3 << " mm\n";
  TVector3 initPos(0, gyroradius, 0);
  TString trackFile{
      "/home/sjones/work/qtnm/Cavities/CavityFieldAmps/track.root"};
  ElectronTrajectoryGen traj(trackFile, field, initPos, eVel, 5e-12, 1e-6, 0,
                             2 * R_E / (3 * TMath::C()));

  // Now reopen the track file and generate some amplitudes
  TFile fTrack(trackFile, "read");
  TTreeReader reader("tree", &fTrack);
  TTreeReaderValue<double> time(reader, "time");
  TTreeReaderValue<double> xPos(reader, "xPos");
  TTreeReaderValue<double> yPos(reader, "yPos");
  TTreeReaderValue<double> zPos(reader, "zPos");
  TTreeReaderValue<double> xVel(reader, "xVel");
  TTreeReaderValue<double> yVel(reader, "yVel");
  TTreeReaderValue<double> zVel(reader, "zVel");

  auto grEn1Real = new TGraph();
  setGraphAttr(grEn1Real);
  auto grEn2Real = new TGraph();
  setGraphAttr(grEn2Real);

  /*
  auto grXPos = new TGraph();
  setGraphAttr(grXPos);
  auto grYPos = new TGraph();
  setGraphAttr(grYPos);
  auto grZPos = new TGraph();
  setGraphAttr(grZPos);
  */
  auto grTTRet = new TGraph();
  setGraphAttr(grTTRet);
  grTTRet->SetTitle("; t [s]; t_{r} [s]");
  TVector3 readoutPos(cavityRadius, 0, 0);  // Hypothetical readout position

  // Loop through entries
  while (reader.Next()) {
    TVector3 r0(*xPos, *yPos, *zPos);
    ComplexVector3 J(*xVel, *yVel, *zVel);
    J *= -TMath::Qe();
    std::complex<double> denom{
        kTE111 * kTE111 -
        k0 * k0 * (1.0 + std::complex<double>(1, -1) / cavityQ)};
    std::complex<double> en1(0, -MU0 * cycFreq * TMath::TwoPi());
    std::complex<double> en2(0, -MU0 * cycFreq * TMath::TwoPi());
    ComplexVector3 eField1{cav.GetModeEField(r0, CircularCavity::kTE,
                                             normalisation1, 1, 1, 1, true, 0)};
    ComplexVector3 eField2{cav.GetModeEField(
        r0, CircularCavity::kTE, normalisation2, 1, 1, 1, false, tTE111 / 4)};
    en1 *= J.Dot(eField1) / denom;
    en2 *= J.Dot(eField2) / denom;
    grEn1Real->SetPoint(grEn1Real->GetN(), *time, en1.real());
    grEn2Real->SetPoint(grEn2Real->GetN(), *time, en2.real());
    grTTRet->SetPoint(grTTRet->GetN(),
                      CalcTimeFromRetardedTime(readoutPos, r0, *time), *time);
  }
  fTrack.Close();

  outPlots.cd();
  grEn1Real->GetXaxis()->SetRangeUser(0, 2e-10);
  grEn2Real->GetXaxis()->SetRangeUser(0, 2e-10);
  grEn2Real->SetLineColor(kRed);
  grEn1Real->Write("grEn1Real");
  grEn2Real->Write("grEn2Real");
  grTTRet->Write("grTTRet");

  // Now make some graphs taking into accout propagation time
  auto spTTRet = new TSpline3("spTTRet", grTTRet);
  auto spEn1 = new TSpline3("spEn1", grEn1Real);
  auto spEn2 = new TSpline3("spEn2", grEn2Real);
  // Equivalent graphs accounting for propagation time
  auto grEn1Real_RO = new TGraph();
  setGraphAttr(grEn1Real_RO);
  grEn1Real_RO->SetTitle("; t [s]; e_{n}");
  auto grEn2Real_RO = new TGraph();
  setGraphAttr(grEn2Real_RO);
  grEn2Real_RO->SetTitle("; t [s]; e_{n}");
  grEn2Real_RO->SetLineColor(kRed);

  // Need to work out the first time that we should take these splines from
  const double timeToStart{grTTRet->GetPointX(0)};
  for (uint i{0}; i < grEn1Real->GetN(); i++) {
    if (grEn1Real->GetPointX(i) < timeToStart) {
      grEn1Real_RO->SetPoint(grEn1Real_RO->GetN(), grEn1Real->GetPointX(i), 0);
      continue;
    } else {
      double tRet{spTTRet->Eval(grEn1Real->GetPointX(i))};
      grEn1Real_RO->SetPoint(grEn1Real_RO->GetN(), grEn1Real->GetPointX(i),
                             spEn1->Eval(tRet));
    }
  }
  delete spEn1;

  for (uint i{0}; i < grEn2Real->GetN(); i++) {
    if (grEn2Real->GetPointX(i) < timeToStart) {
      grEn2Real_RO->SetPoint(grEn2Real_RO->GetN(), grEn2Real->GetPointX(i), 0);
      continue;
    } else {
      double tRet{spTTRet->Eval(grEn2Real->GetPointX(i))};
      grEn2Real_RO->SetPoint(grEn2Real_RO->GetN(), grEn2Real->GetPointX(i),
                             spEn2->Eval(tRet));
    }
  }
  delete spEn2;

  outPlots.cd();
  grEn1Real_RO->Write("grEn1Real_RO");
  grEn2Real_RO->Write("grEn2Real_RO");
  delete spTTRet;

  // Now generate power spectra
  auto grEn1PowerSpec = MakePowerSpectrumPeriodogram(grEn1Real_RO);
  auto grEn2PowerSpec = MakePowerSpectrumPeriodogram(grEn2Real_RO);
  grEn1PowerSpec->SetLineColor(kRed);
  outPlots.cd();
  grEn1PowerSpec->Write("grEn1PowerSpec");
  grEn2PowerSpec->Write("grEn2PowerSpec");

  outPlots.Close();
  return 0;
}