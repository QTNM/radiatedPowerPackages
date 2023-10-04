/// WR42Jackson.cxx

#include <complex>
#include <iostream>
#include <memory>

#include "BasicFunctions/BasicFunctions.h"
#include "BasicFunctions/ComplexVector3.h"
#include "BasicFunctions/Constants.h"
#include "ElectronDynamics/QTNMFields.h"
#include "ElectronDynamics/TrajectoryGen.h"
#include "TF1.h"
#include "TFile.h"
#include "TMath.h"
#include "TString.h"
#include "TTree.h"
#include "Waveguides/RectangularWaveguide.h"

using namespace rad;

double GetEModeNormFactor(RectangularWaveguide *wv, int m, int n, double freq,
                          int nSurfPnts) {
  double area{(wv->GetLongDimension() / double(nSurfPnts)) *
              (wv->GetShortDimension() / double(nSurfPnts))};
  double sum{0.0};

  for (int ix = 0; ix < nSurfPnts; ix++) {
    double thisx{-wv->GetLongDimension() / 2.0 +
                 wv->GetLongDimension() / (2.0 * double(nSurfPnts)) +
                 wv->GetLongDimension() * double(ix) / double(nSurfPnts)};
    for (int iy = 0; iy < nSurfPnts; iy++) {
      double thisy{-wv->GetShortDimension() / 2.0 +
                   wv->GetShortDimension() / (2.0 * double(nSurfPnts)) +
                   wv->GetShortDimension() * double(iy) / double(nSurfPnts)};
      TVector3 surfacePos{thisx, thisy, 0.0};
      ComplexVector3 eTransReal{wv->GetNormalisedEField(
          RectangularWaveguide::kTE, m, n, surfacePos, freq * 2 * TMath::Pi())};
      eTransReal.SetZ(0.0);
      sum += (eTransReal.Dot(eTransReal)).real() * area;
    }  // Loop over y points
  }    // Loop over x points

  return sum;
}

double GetHModeNormFactor(RectangularWaveguide *wv, int m, int n, double freq,
                          int nSurfPnts) {
  double area{(wv->GetLongDimension() / double(nSurfPnts)) *
              (wv->GetShortDimension() / double(nSurfPnts))};
  double sum{0.0};

  for (int ix = 0; ix < nSurfPnts; ix++) {
    double thisx{-wv->GetLongDimension() / 2.0 +
                 wv->GetLongDimension() / (2.0 * double(nSurfPnts)) +
                 wv->GetLongDimension() * double(ix) / double(nSurfPnts)};
    for (int iy = 0; iy < nSurfPnts; iy++) {
      double thisy{-wv->GetShortDimension() / 2.0 +
                   wv->GetShortDimension() / (2.0 * double(nSurfPnts)) +
                   wv->GetShortDimension() * double(iy) / double(nSurfPnts)};
      TVector3 surfacePos{thisx, thisy, 0.0};
      ComplexVector3 hTransReal{wv->GetNormalisedHField(
          RectangularWaveguide::kTE, m, n, surfacePos, freq * 2 * TMath::Pi())};
      hTransReal.SetZ(0.0);
      sum += (hTransReal.Dot(hTransReal)).real() * area;
    }  // Loop over y points
  }    // Loop over x points

  return sum;
}

//////////////////////////// Field amplitudes
////////////////////////////////////////////
std::complex<double> GetPositiveAmp(RectangularWaveguide *wv, int m, int n,
                                    double freq, TVector3 pos, TVector3 vel) {
  double waveImp{wv->GetModeImpedance(RectangularWaveguide::kTE, m, n,
                                      freq * 2 * TMath::Pi())};
  TVector3 j{-TMath::Qe() * vel};
  ComplexVector3 jComplex{j};

  ComplexVector3 eTrans{wv->GetNormalisedEField(RectangularWaveguide::kTE, m, n,
                                                pos, freq * 2 * TMath::Pi())};
  eTrans.SetZ(std::complex<double>{0.0, 0.0});
  ComplexVector3 eAxial{wv->GetNormalisedEField(RectangularWaveguide::kTE, m, n,
                                                pos, freq * 2 * TMath::Pi())};
  eAxial.SetX(std::complex<double>{0.0, 0.0});
  eAxial.SetY(std::complex<double>{0.0, 0.0});
  ComplexVector3 subVec{eTrans - eAxial};

  std::complex<double> APlus = (subVec).Dot(jComplex) * (-waveImp / 2.0);
  return APlus;
}

std::complex<double> GetNegativeAmp(RectangularWaveguide *wv, int m, int n,
                                    double freq, TVector3 pos, TVector3 vel) {
  double waveImp{wv->GetModeImpedance(RectangularWaveguide::kTE, m, n,
                                      freq * 2 * TMath::Pi())};
  TVector3 j{-TMath::Qe() * vel};
  ComplexVector3 jComplex{j};

  ComplexVector3 eTrans{wv->GetNormalisedEField(RectangularWaveguide::kTE, m, n,
                                                pos, freq * 2 * TMath::Pi())};
  eTrans.SetZ(std::complex<double>{0.0, 0.0});
  ComplexVector3 eAxial{wv->GetNormalisedEField(RectangularWaveguide::kTE, m, n,
                                                pos, freq * 2 * TMath::Pi())};
  eAxial.SetX(std::complex<double>{0.0, 0.0});
  eAxial.SetY(std::complex<double>{0.0, 0.0});
  ComplexVector3 subVec{eTrans + eAxial};

  std::complex<double> AMinus = (subVec).Dot(jComplex) * (-waveImp / 2.0);
  return AMinus;
}

double CalculatePowerPlus(RectangularWaveguide *wv, int m, int n, double freq,
                          TVector3 ePos, TVector3 eVel, double integralZPos,
                          int nSurfPnts) {
  // Firstly get the field amplitudes
  std::complex<double> aPlus{GetPositiveAmp(wv, m, n, freq, ePos, eVel)};
  double totalPowerPlus{0.0};
  double area{(wv->GetLongDimension() / double(nSurfPnts)) *
              (wv->GetShortDimension() / double(nSurfPnts))};

  for (int ix = 0; ix < nSurfPnts; ix++) {
    double thisx{-wv->GetLongDimension() / 2.0 +
                 wv->GetLongDimension() / (2.0 * double(nSurfPnts)) +
                 wv->GetLongDimension() * double(ix) / double(nSurfPnts)};
    for (int iy = 0; iy < nSurfPnts; iy++) {
      double thisy{-wv->GetShortDimension() / 2.0 +
                   wv->GetShortDimension() / (2.0 * double(nSurfPnts)) +
                   wv->GetShortDimension() * double(iy) / double(nSurfPnts)};
      TVector3 surfacePos{thisx, thisy, integralZPos};

      ComplexVector3 ETotal{wv->GetNormalisedEField(
          RectangularWaveguide::kTE, m, n, surfacePos, freq * 2 * TMath::Pi())};
      ETotal *= aPlus;
      ComplexVector3 HTotal{wv->GetNormalisedHField(
          RectangularWaveguide::kTE, m, n, surfacePos, freq * 2 * TMath::Pi())};
      HTotal *= aPlus;

      totalPowerPlus += (ETotal.Cross(HTotal.Conj())).Z().real() * area;
    }
  }

  return totalPowerPlus;
}

int main(int argc, char *argv[]) {
  TString outputFile{argv[1]};
  auto fout = std::make_unique<TFile>(outputFile, "RECREATE");

  // Generate the field
  // This is chosen so we get the correct frequency
  double rCoil{3e-3};
  double trapDepth{4e-3};
  UniformField *field = new UniformField(0.964);

  TString trackFile{"/home/sjones/work/qtnm/outputs/WR42Jackson/track.root"};
  // Calculate the central frequency
  const double centralBField{
      field->evaluate_field_at_point(TVector3(0, 0, 0)).Mag()};
  const double centralFreq{26e9};  // GHz
  std::cout << "Central frequency = " << centralFreq / 1e9 << " GHz\n"
            << std::endl;

  const double ke{
      (TMath::Qe() * centralBField / (2 * TMath::Pi() * centralFreq) - ME) *
      pow(TMath::C(), 2) / TMath::Qe()};
  std::cout << "KE = " << ke << std::endl;
  const double speed{GetSpeedFromKE(ke, ME)};
  const TVector3 velocity{speed, 0, 0};
  const double gyroradius{GetGyroradius(
      velocity, field->evaluate_field_at_point(TVector3(0, 0, 0)), ME)};
  std::cout << "Gyroradius = " << gyroradius * 1e3 << " mm" << std::endl;

  // Project 8 WR42 waveguide used the TE10 mode
  const double WR42Side1{10.668e-3};
  const double WR42Side2{4.318e-3};
  const double waveguideLength{0.05};  // Don't think this is very important
  TVector3 probePosition(0, WR42Side2 * 0.2, 0);
  RectangularWaveguide *WR42 = new RectangularWaveguide(
      WR42Side1, WR42Side2, waveguideLength, probePosition);

  // Check fields are normalised
  const int nSurfPnts{30};
  double integralE{GetEModeNormFactor(WR42, 1, 0, centralFreq, nSurfPnts)};
  double integralH{GetHModeNormFactor(WR42, 1, 0, centralFreq, nSurfPnts)};
  double waveImp{WR42->GetModeImpedance(RectangularWaveguide::kTE, 1, 0,
                                        centralFreq * 2 * TMath::Pi())};
  std::cout << "Integral E = " << integralE << std::endl;  // Should be 1
  std::cout << "Integral H = " << integralH
            << std::endl;  // Should equal to the below
  std::cout << "1 over wave impedance squared = " << 1.0 / (waveImp * waveImp)
            << std::endl;

  std::cout << "\n";
  const int nOffsets{51};
  const double maxOffset{4.8e-3};   // metres
  const double minOffset{-4.8e-3};  // metres

  TGraph *grCollectedPower = new TGraph();
  setGraphAttr(grCollectedPower);
  grCollectedPower->SetMarkerStyle(20);
  grCollectedPower->SetTitle(
      "Collected power (unscaled): WR42 waveguide; X [mm]; Collected power "
      "[fW]");
  TGraph *grCollectedPowerCrude = new TGraph();
  setGraphAttr(grCollectedPowerCrude);
  grCollectedPowerCrude->SetMarkerStyle(20);
  grCollectedPowerCrude->SetTitle(
      "Collected power (unscaled): WR42 waveguide; X [mm]; Collected power "
      "[fW]");
  TGraph *grCollectedFrac = new TGraph();
  setGraphAttr(grCollectedFrac);
  grCollectedFrac->SetMarkerStyle(20);
  grCollectedFrac->SetTitle(
      "Collected power fraction: WR42 waveguide; X [mm]; Fraction of total "
      "power");

  for (int iOffset = 0; iOffset < nOffsets; iOffset++) {
    double testOffset{minOffset + double(iOffset) * (maxOffset - minOffset) /
                                      double(nOffsets - 1)};
    std::cout << testOffset * 1e3 << " mm" << std::endl;

    // Generate the electron trajectory
    ElectronTrajectoryGen traj(trackFile, field,
                               TVector3(testOffset, -gyroradius, 0.0),
                               TVector3(speed, 0, 0), 1e-12, 1e-9);

    // Now read in the file again and calculate some mode amplitudes as a
    // function of time
    auto fin = std::make_unique<TFile>(trackFile, "READ");
    TTree *tree = (TTree *)fin->Get("tree");
    double time{};
    double xPos{}, yPos{}, zPos{};
    double xVel{}, yVel{}, zVel{};
    tree->SetBranchAddress("time", &time);
    tree->SetBranchAddress("xPos", &xPos);
    tree->SetBranchAddress("yPos", &yPos);
    tree->SetBranchAddress("zPos", &zPos);
    tree->SetBranchAddress("xVel", &xVel);
    tree->SetBranchAddress("yVel", &yVel);
    tree->SetBranchAddress("zVel", &zVel);

    const double volumeLength{0.005};

    TGraph *grPowerPlus = new TGraph();
    setGraphAttr(grPowerPlus);
    grPowerPlus->SetTitle(
        Form("x = %.2f: P^{+}; Time [s]; P^{+}", testOffset * 1e3));
    grPowerPlus->SetMarkerStyle(20);
    TGraph *grCrudePlus = new TGraph();
    setGraphAttr(grCrudePlus);
    grCrudePlus->SetTitle(
        Form("x = %.2f: P^{+}; Time [s]; P^{+}", testOffset * 1e3));
    TGraph *grAPlus = new TGraph();
    setGraphAttr(grAPlus);
    grAPlus->SetTitle(
        Form("x = %.2f: A^{+}; Time [s]; A^{+}", testOffset * 1e3));

    TGraph *grEnergy = new TGraph();
    setGraphAttr(grEnergy);
    grEnergy->SetTitle("Electron energy; Time [s]; Kinetic energy [J]");

    // Loop over entries
    for (int iE = 0; iE < tree->GetEntries(); iE++) {
      tree->GetEntry(iE);
      TVector3 thePos{xPos, yPos, zPos};
      TVector3 theVel{xVel, yVel, zVel};

      double beta{theVel.Mag() / TMath::C()};
      double gamma{1.0 / sqrt(1 - beta * beta)};
      double thiske{(gamma - 1.0) * ME * TMath::C() * TMath::C()};
      grEnergy->SetPoint(iE, time, thiske);

      std::complex<double> ampPlus{
          GetPositiveAmp(WR42, 1, 0, centralFreq, thePos, theVel)};
      double totalPowerPlus{CalculatePowerPlus(
          WR42, 1, 0, centralFreq, thePos, theVel, volumeLength, nSurfPnts)};
      grAPlus->SetPoint(grAPlus->GetN(), time, ampPlus.real());
      grCrudePlus->SetPoint(grCrudePlus->GetN(), time,
                            pow(abs(ampPlus), 2) * 1e15 / waveImp);
      grPowerPlus->SetPoint(grPowerPlus->GetN(), time, totalPowerPlus * 1e15);
    }  // Loop over electron trajectory entries
    fin->Close();

    fout->cd();
    grAPlus->Write(Form("grAPlus%d", iOffset));
    grCrudePlus->Write(Form("grCrudePlus%d", iOffset));
    grPowerPlus->Write(Form("grPowerPlus%d", iOffset));

    // Determine the radiated power
    TF1 *fEnergyLoss = new TF1("fEnergyLoss", "pol1", 0, 1e-9);
    grEnergy->Fit(fEnergyLoss);
    double radiatedPower{-1 * fEnergyLoss->GetParameter(1)};
    std::cout << "Radiated power is " << radiatedPower * 1e15 << " fW"
              << std::endl;

    double timeAvgPower{0};
    double timeAvgPowerCrude{0};
    for (int n = 0; n < grCrudePlus->GetN(); n++) {
      timeAvgPower += grPowerPlus->GetPointY(n);
      timeAvgPowerCrude += grCrudePlus->GetPointY(n);
    }
    timeAvgPower /= double(grPowerPlus->GetN());
    timeAvgPowerCrude /= double(grCrudePlus->GetN());
    grCollectedPower->SetPoint(grCollectedPower->GetN(), testOffset * 1e3,
                               timeAvgPower);
    grCollectedPowerCrude->SetPoint(grCollectedPowerCrude->GetN(),
                                    testOffset * 1e3, timeAvgPowerCrude);
    grCollectedFrac->SetPoint(grCollectedFrac->GetN(), testOffset * 1e3,
                              timeAvgPower / (radiatedPower * 1e15));

    delete fEnergyLoss;
    delete grEnergy;
  }  // Loop over offsets

  fout->cd();
  grCollectedPower->Write("grCollectedPower");
  grCollectedPowerCrude->Write("grCollectedPowerCrude");
  grCollectedFrac->Write("grCollectedFrac");

  fout->Close();
  return 0;
}
