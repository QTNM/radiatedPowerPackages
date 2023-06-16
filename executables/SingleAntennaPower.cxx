/*
  SingleAntennaPower.cxx

  Measurement of the power collected by a single dipole antenna
*/

#include "Antennas/HalfWaveDipole.h"
#include "BasicFunctions/BasicFunctions.h"
#include "BasicFunctions/Constants.h"
#include "ElectronDynamics/QTNMFields.h"
#include "ElectronDynamics/TrajectoryGen.h"
#include "FieldClasses/FieldClasses.h"
#include "FieldClasses/FieldPointNR.h"
#include "SignalProcessing/InducedVoltage.h"
#include "TFile.h"
#include "TGraph.h"
#include "TMath.h"
#include "TString.h"
#include "TVector3.h"

using namespace rad;
using std::cout;
using std::endl;

double LarmorPowerNR(double a) {
  return MU0 * pow(TMath::Qe() * a, 2) / (6 * TMath::Pi() * TMath::C());
}

int main(int argc, char *argv[]) {
  TString outputFile{argv[1]};
  TFile *fout = new TFile(outputFile, "recreate");

  // Field details
  const double fieldMag{1.0};  // Tesla
  UniformField *field = new UniformField(fieldMag);
  const TVector3 centralField{
      field->evaluate_field_at_point(TVector3(0, 0, 0))};

  // Electron details
  const double electronKE{18.6e3};  // eV
  const double electronSpeed{GetSpeedFromKE(electronKE, ME)};
  const double tau{2 * R_E / (3 * TMath::C())};
  TVector3 vel0(electronSpeed, 0, 0);
  const double gyroradius{GetGyroradius(vel0, centralField, ME)};  // metres
  TVector3 pos0(0, -gyroradius, 0);
  const double centralFreq{CalcCyclotronFreq(electronKE, centralField.Mag())};
  const double acc0{electronSpeed * 2 * TMath::Pi() * centralFreq};
  const double gamma{1 / sqrt(1 - pow(electronSpeed / TMath::C(), 2))};

  const double radiatedPower{
      MU0 *
      pow(TMath::Qe() * 2 * TMath::Pi() * centralFreq * electronSpeed, 2) *
      pow(gamma, 4) / (6 * TMath::Pi() * TMath::C())};
  const double radiatedPowerNR{LarmorPowerNR(acc0)};
  cout << "Radiated power (R, NR) = " << radiatedPower * 1e15 << " fW,\t"
       << radiatedPowerNR << " fW\n";

  // Trajectory generating details
  const double simTime{1e-6};       // seconds
  const double simStepSize{1e-12};  // seconds
  TString trackFile{
      "/home/sjones/work/qtnm/outputs/SingleAntennaPower/track.root"};
  ElectronTrajectoryGen traj(trackFile, field, pos0, vel0, simStepSize, simTime,
                             0.0, tau);

  const double xAxisDeg{1500};
  const double xAxisCycles{xAxisDeg / 360};
  const double xAxisTime{xAxisCycles * (1 / centralFreq)};

  // Now we've generated the trajectory we want to measure the Poynting vector
  // at an antenna - use a dipole for simplicity
  // Firstly try an antenna 2cm away from the electron
  const double antennaRadius{0.02};  // metres
  TVector3 antPos(antennaRadius, 0, 0);
  TVector3 antXDir(1, 0, 0);
  TVector3 antZDir(0, 1, 0);
  HalfWaveDipole *antenna =
      new HalfWaveDipole(antPos, antXDir, antZDir, centralFreq);
  FieldPoint fp(trackFile, antenna);
  fp.GenerateFields(0, simTime);
  TGraph *grS = fp.GetPoyntingMagTimeDomain(true);
  // Calculate the time average of the Poynting vector
  cout << "\n";
  double avgS{0};
  for (int n{0}; n < grS->GetN(); n++) {
    avgS += grS->GetPointY(n);
  }
  avgS /= double(grS->GetN());
  cout << "Time averaged Poynting vector = " << avgS << " W m^-2\n";

  fout->cd();
  grS->GetXaxis()->SetRangeUser(1e-9, 1e-9 + xAxisTime);
  grS->Write("grS");

  // Now get the collected power as a function of time (using the effective
  // area)
  cout << "\n";
  TGraph *grPowerAEff = fp.GetAntennaPowerTimeDomain(true);
  // Calculate time averaged power
  double avgP{0};
  for (int n{0}; n < grPowerAEff->GetN(); n++) {
    avgP += grPowerAEff->GetPointY(n);
  }
  avgP /= double(grPowerAEff->GetN());
  cout << "Time averaged power (effective area) = " << avgP * 1e15 << " fW\n";
  cout << "Efficiency = " << avgP * 100 / radiatedPower << "%\n";

  fout->cd();
  grPowerAEff->GetYaxis()->SetRangeUser(0, 0.03e-15);
  grPowerAEff->GetXaxis()->SetRangeUser(1e-9, 1e-9 + xAxisTime);
  grPowerAEff->Write("grPowerAEff");

  // Try this for the non-relativistic case also
  cout << "\n";
  FieldPointNR fpNR(trackFile, antenna);
  fpNR.GenerateFields(0, simTime);
  TGraph *grPowerAEffNR{fpNR.GetAntennaPowerTimeDomain(true)};

  double avgPNR{0};
  for (int n{0}; n < grPowerAEffNR->GetN(); n++) {
    avgPNR += grPowerAEffNR->GetPointY(n);
  }
  avgPNR /= double(grPowerAEffNR->GetN());
  cout << "Time averaged non-relativistic power (effective area) = "
       << avgPNR * 1e15 << " fW\n";
  cout << "Efficiency = " << avgPNR * 100 / radiatedPowerNR << "%\n";
  fout->cd();
  grPowerAEffNR->GetXaxis()->SetRangeUser(1e-9, 1e-9 + xAxisTime);
  grPowerAEffNR->GetYaxis()->SetRangeUser(0, 8e-18);
  grPowerAEffNR->Write("grPowerAEffNR");

  // Now get the induced voltage and see if we get a similar result
  cout << "\n";
  InducedVoltage iv(trackFile, antenna, true);
  iv.GenerateVoltage(0, simTime);
  TGraph *grV = iv.GetVoltageGraph();

  fout->cd();
  grV->Write("grV");

  TGraph *grPowerV = new TGraph();
  setGraphAttr(grPowerV);
  grPowerV->SetMarkerStyle(20);
  grPowerV->SetTitle("Single dipole; Time [s]; Power [W]");
  double avgPV{0};
  for (int n{0}; n < grV->GetN(); n++) {
    double power{grV->GetPointY(n) * grV->GetPointY(n) / 73};
    avgPV += power;
    grPowerV->SetPoint(n, grV->GetPointX(n), power);
  }
  avgPV /= double(grPowerV->GetN());
  cout << "Time averaged power (voltage) = " << avgPV * 1e15 << " fW\n";

  fout->cd();
  grPowerV->GetXaxis()->SetRangeUser(1e-9, 1e-9 + xAxisTime);
  grPowerV->GetYaxis()->SetRangeUser(0, 0.03e-15);
  grPowerV->Write("grPowerV");

  // Make power spectrum for this voltage
  TGraph *grVPgram = iv.GetPowerPeriodogram(73);
  fout->cd();
  grVPgram->Write("grVPgram");
  double pgramSum{0};
  double pgramSumFRange{0};  // Only count frequency range we'd detect
  for (int n{0}; n < grVPgram->GetN(); n++) {
    pgramSum += grVPgram->GetPointY(n);
    if (grVPgram->GetPointX(n) > centralFreq - 500e6 &&
        grVPgram->GetPointX(n) < centralFreq + 500e6) {
      pgramSumFRange += grVPgram->GetPointY(n);
    }
  }
  cout << "Periodogram power sum (in range of interest)= " << pgramSum * 1e15
       << " (" << pgramSumFRange * 1e15 << ") fW\n";
  cout << "Efficiency = " << pgramSum / radiatedPower * 100 << "%\n";

  fout->Close();
  delete fout;
  return 0;
}
