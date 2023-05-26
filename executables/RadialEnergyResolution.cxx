/*
  RadialEnergyResolution.cxx

  Test how limited knowledge of radial magnetic field affects energy resolution
*/

#include "BasicFunctions/Constants.h"
#include "BasicFunctions/BasicFunctions.h"

#include "ElectronDynamics/QTNMFields.h"
#include "ElectronDynamics/ComsolFields.h"

#include "TFile.h"
#include "TGraph.h"
#include "TH1.h"
#include "TMath.h"
#include "TRandom3.h"
#include "TString.h"
#include "TVector3.h"
#include "THStack.h"
#include "TF1.h"

#include <iostream>
#include <memory>

using namespace rad;
using std::unique_ptr;

double SkewedGaussian_f(double *x, double *par)
{
  // return

  double xx{x[0]};
  double A{par[0]};
  double mu{par[1]};
  double sigma{par[2]};
  double alpha{par[3]};
  return SkewedGaussian(xx, A, mu, sigma, alpha);
}

int main(int argc, char *argv[])
{
  TString outputFile{argv[1]};
  TFile *fout{new TFile(outputFile, "recreate")};

  // Field specifications
  const double desiredB{1}; // Tesla
  const double radius{0.2}; // metres
  const double length{1.5}; // metres
  const double current{20}; // Amps
  const double n{desiredB / (MU0 * current)};
  unique_ptr<SolenoidField> field =
      std::make_unique<SolenoidField>(radius, length, current, n);
  const double centralField{field->evaluate_field_magnitude(TVector3(0, 0, 0))};

  // Now grab one of the COMSOL fields to have a look at that
  unique_ptr<ComsolField> cfield =
      std::make_unique<ComsolField>("/home/sjones/work/qtnm/radiatedPowerPackages/files/fieldmap_50A_170mmx1000mm.csv");
  const double centralFieldComsol{cfield->evaluate_field_magnitude(TVector3(0, 0, 0))};

  // First plot just the radial field and the variation in it
  TGraph *grDeltaB{new TGraph()};
  setGraphAttr(grDeltaB);
  grDeltaB->SetTitle(Form("B_{central} = %.1f T", desiredB));
  grDeltaB->GetXaxis()->SetTitle("R [m]");
  grDeltaB->GetYaxis()->SetTitle("#Delta B / B");

  TGraph *grDeltaBComsol{new TGraph()};
  setGraphAttr(grDeltaBComsol);
  grDeltaBComsol->SetTitle("COMSOL field");
  grDeltaBComsol->GetXaxis()->SetTitle("R [m]");
  grDeltaBComsol->GetYaxis()->SetTitle("#Delta B / B");

  // Some kinematic info that we need
  const double eKE{18.5e3}; // eV
  const double gamma{(ME * TMath::C() * TMath::C() + eKE * TMath::Qe()) / (ME * TMath::C() * TMath::C())};
  std::cout << "Gamma = " << gamma << std::endl;

  const double minRadius{0.001};
  const double maxRadius{0.085};
  const int nRPnts{237};

  for (int iR{0}; iR < nRPnts; iR++)
  {
    double testRadius{minRadius + (maxRadius - minRadius) * double(iR) / double(nRPnts - 1)};
    double testField{field->evaluate_field_magnitude(TVector3(0, testRadius, 0))};
    grDeltaB->SetPoint(iR, testRadius, (testField - centralField) / centralField);
    double testFieldComsol{cfield->evaluate_field_magnitude(TVector3(0.001, testRadius, 0.001))};
    grDeltaBComsol->SetPoint(iR, testRadius,
                             (testFieldComsol - centralFieldComsol) / centralFieldComsol);
  }

  fout->cd();
  grDeltaB->Write("grDeltaB");
  grDeltaBComsol->Write("grDeltaBComsol");

  // Randomly distribute points (uniform density) and check what the field is
  TH1D *hB{new TH1D("hB", "B field; B_{init} [T]; N_{electrons}",
                    100, centralField, centralField * (1 + 0.9e-3))};
  SetHistAttr(hB);
  TRandom3 *rand{new TRandom3()};
  const int nElectrons{40000};
  for (int iE{0}; iE < nElectrons; iE++)
  {
    double eR{0.1 * sqrt(rand->Uniform(0, 1))};
    double bField{field->evaluate_field_magnitude(TVector3(eR, 0, 0))};
    hB->Fill(bField);
  }

  // Assume we our field perfectly matches our model
  // Only uncertainty is on electron position
  // Test what kind of B field resolutions we get
  const double posResMin{0.001};
  const double posResMax{0.05};
  const int nResPnts{30};
  double logRDiff{(log10(posResMax) - log10(posResMin)) / double(nResPnts - 1)};
  THStack *hs = new THStack("hs", "B field resolution as a function of radial position resolution; #Delta B / B; N_{e}");

  TGraph *grResRResB{new TGraph()};
  setGraphAttr(grResRResB);
  grResRResB->SetMarkerStyle(20);
  grResRResB->SetTitle("B field resolution from radial position resolution");
  grResRResB->GetXaxis()->SetTitle("#sigma_{R} [mm]");
  grResRResB->GetYaxis()->SetTitle("#sigma_{B} / B");

  TGraph *grResRResE{new TGraph()};
  setGraphAttr(grResRResE);
  grResRResE->SetMarkerStyle(20);
  grResRResE->SetTitle("Energy resolution from radial position resolution");
  grResRResE->GetXaxis()->SetTitle("#sigma_{R} [mm]");
  grResRResE->GetYaxis()->SetTitle("#sigma_{E} / E");

  for (int n{0}; n < nResPnts; n++)
  {
    double res{posResMin * pow(10, double(n) * logRDiff)};
    TH1D *hBRes = new TH1D(Form("hBRes%d", n), Form("#sigma_{R} = %.4f", res), 400, -0.0015, 0.0015);
    SetHistAttr(hBRes);
    hBRes->GetXaxis()->SetTitle("#Delta B / B");
    hBRes->GetYaxis()->SetTitle("N_{e}");
    hBRes->SetLineColor(52 + n * 2);

    TF1 *fSGaus = new TF1(Form("fSGaus%d", n), SkewedGaussian_f, -0.002, 0.002, 4);
    fSGaus->SetParName(0, "Scale");
    fSGaus->SetParName(1, "#mu");
    fSGaus->SetParName(2, "#sigma");
    fSGaus->SetParName(3, "Skewness");

    for (int iE{0}; iE < nElectrons; iE++)
    {
      double rTrue{0.1 * sqrt(rand->Uniform(0, 1))};
      double bTrue{field->evaluate_field_magnitude(TVector3(rTrue, 0, 0))};
      double rReco{rand->Gaus(rTrue, res)};
      double bReco{field->evaluate_field_magnitude(TVector3(rReco, 0, 0))};
      hBRes->Fill((bReco - bTrue) / bTrue);
    }
    hs->Add(hBRes);

    fSGaus->SetParameter(0, hBRes->GetMaximum());
    fSGaus->SetParameter(1, hBRes->GetBinCenter(hBRes->GetMaximumBin()));
    fSGaus->SetParameter(2, hBRes->GetRMS());
    fSGaus->SetParameter(3, 0);
    fSGaus->SetLineColor(kBlack);

    TFitResultPtr fresult = hBRes->Fit(fSGaus, "r");
    int fitStatus = fresult;
    std::cout << "Fit result = " << fitStatus << std::endl;

    if (fitStatus == 0)
    {
      grResRResB->SetPoint(n, res * 1e3, fSGaus->GetParameter(2));
      grResRResE->SetPoint(n, res * 1e3,
                           fSGaus->GetParameter(2) * gamma / (gamma - 1));
    }

    fout->cd();
    hBRes->Write();
    fSGaus->Write();
    delete fSGaus;
    std::cout << "\n";
  }

  fout->cd();
  hs->Write();
  hB->Write();
  grResRResB->Write("grResRResB");
  grResRResE->Write("grResRResE");

  fout->Close();
  delete fout;
  return 0;
}