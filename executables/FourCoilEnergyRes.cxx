/*
  FourCoilEnergyRes.cxx

  Generate tracks in multiple different configurations
*/

#include <iostream>

#include "BasicFunctions/BasicFunctions.h"
#include "BasicFunctions/Constants.h"
#include "ElectronDynamics/QTNMFields.h"
#include "TF1.h"
#include "TFile.h"
#include "TGraph.h"
#include "TString.h"

using namespace rad;

double UniformBkg(double *x, double *par) {
  // Trap coil definitions
  const double rCoil{15e-3};                                            // m
  const double outerCoilOffset{25e-3};                                  // m
  const double outerCoilTrapDepth{5.5e-3};                              // T
  const double outerCoilCurrent{2 * outerCoilTrapDepth * rCoil / MU0};  // A
  const double innerCoilOffset{19e-3};                                  // m
  const double innerCoilTrapDepth{3.1e-3};                              // T
  const double innerCoilCurrent{2 * innerCoilTrapDepth * rCoil / MU0};  // A

  // Set up uniform background field first
  const double bkgField{0.7};  // Tesla
  FourCoilBathtub uniformBkg(rCoil, outerCoilCurrent, outerCoilOffset, rCoil,
                             -innerCoilCurrent, innerCoilOffset, bkgField);

  double zPos{x[0] / 1e3};
  TVector3 pos(0, 0, zPos);
  return uniformBkg.evaluate_field_magnitude(pos);
}

double HelmholtzBkg(double *x, double *par) {
  // Trap coil definitions
  const double rCoil{15e-3};                                            // m
  const double outerCoilOffset{25e-3};                                  // m
  const double outerCoilTrapDepth{5.5e-3};                              // T
  const double outerCoilCurrent{2 * outerCoilTrapDepth * rCoil / MU0};  // A
  const double innerCoilOffset{19e-3};                                  // m
  const double innerCoilTrapDepth{3.1e-3};                              // T
  const double innerCoilCurrent{2 * innerCoilTrapDepth * rCoil / MU0};  // A

  // Helmholtz definitions
  const double n{430};
  const double lSolenoid{5e-2};
  const double rSolenoid{200e-3};
  const double bField{0.4912};
  const double iSolenoid{2 * bField *
                         sqrt(pow(lSolenoid / 2, 2) + rSolenoid * rSolenoid) /
                         (MU0 * n * lSolenoid)};  // Amps

  HelmholtzField bkg(rSolenoid, lSolenoid, iSolenoid, n);
  FourCoilField trap(rCoil, outerCoilCurrent, outerCoilOffset, rCoil,
                     -innerCoilCurrent, innerCoilOffset);
  FourCoilHelmholtz combinedField(trap, bkg);

  double zPos{x[0] / 1e3};
  TVector3 pos(0, 0, zPos);
  return combinedField.evaluate_field_magnitude(pos);
}

int main(int argc, char *argv[]) {
  TString outputFile{
      "/home/sjones/work/qtnm/MagnetDesigns/EnergyResolution/totalFields.root"};
  TFile fout(outputFile, "recreate");

  TF1 fUniformBkg("fUniformBkg", UniformBkg, -80, 80, 0);
  fUniformBkg.SetNpx(1000);
  fUniformBkg.SetLineWidth(3);
  fUniformBkg.SetTitle("Uniform background; z [mm]; |B| [T]");

  TF1 fHelmholtzBkg("fHelmholtzBkg", HelmholtzBkg, -80, 80, 0);
  fHelmholtzBkg.SetNpx(1000);
  fHelmholtzBkg.SetLineWidth(3);
  fHelmholtzBkg.SetTitle("Helmholtz background; z [mm]; |B| [T]");
  fHelmholtzBkg.SetLineColor(kBlue);

  fout.cd();
  fUniformBkg.Write();
  fHelmholtzBkg.Write();

  fout.Close();
  return 0;
}