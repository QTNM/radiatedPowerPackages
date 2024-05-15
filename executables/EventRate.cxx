/*
  EventRate.cxx

  How many tritium events do we expect to see?
*/

#include <getopt.h>

#include <cmath>
#include <iostream>
#include <string>

#include "BasicFunctions/BasicFunctions.h"
#include "BasicFunctions/Constants.h"
#include "BasicFunctions/TritiumSpectrum.h"
#include "TF1.h"
#include "TFile.h"
#include "TMath.h"
#include "TString.h"

using namespace rad;

/// @brief
/// @param f Frequency in Hz
/// @param B Magnetic field in T
/// @return Electric kinetic energy in eV
double CalcEFromFreq(double f, double B) {
  double eTotJ{TMath::Qe() * B * pow(TMath::C(), 2) / (f * TMath::TwoPi())};
  return (eTotJ - ME * TMath::C() * TMath::C()) / TMath::Qe();
}

double TritiumESpectrum(double *x, double *par) {
  double E{x[0]};
  return TritiumDecayRate(E, 0, 0, 0);
}

double TritiumFSpectrum(double *x, double *par) {
  double f{x[0] * 1e9};
  double E{CalcEFromFreq(f, 0.7)};
  return TritiumDecayRate(E, 0, 0, 0);
}

int main(int argc, char *argv[]) {
  double bField{0.7};    // Tesla
  double deltaF{100e6};  // Hertz
  std::string outputStr{};
  int opt{};
  while ((opt = getopt(argc, argv, ":o:b:f:")) != -1) {
    switch (opt) {
      case 'o':
        outputStr = optarg;
        break;

      case 'b':
        bField = std::stod(optarg);
        break;

      case 'f':
        deltaF = std::stod(optarg);
        break;

      case ':':
        std::cout << "Option needs a value\n";
        break;

      case '?':
        std::cout << "Unknown option: " << optopt << std::endl;
        break;
    }
  }

  TString outputFile{"/home/sjones/work/qtnm/storageSize/plots.root"};
  TFile fout(outputFile, "recreate");

  const double energyBeyondEndpoint{5};
  const double eMax{18.6e3 + energyBeyondEndpoint};    // eV
  const double fMin{CalcCyclotronFreq(eMax, bField)};  // Hz
  const double fMax{fMin + deltaF};                    // Hz
  const double eMin{CalcEFromFreq(fMax, bField)};      // eV
  std::cout << "Emin, Emax = " << eMin << " eV, " << eMax
            << " eV\tDelta E = " << (eMax - eMin) << " eV\n";

  TF1 fESpec("fESpec", TritiumESpectrum, eMin, eMax);
  fESpec.SetTitle("; E [keV]; d#Gamma/dE [s^{-1} eV^{-1}]");
  fESpec.SetLineWidth(3);
  fESpec.SetNpx(500);
  TF1 fFSpec("fFSpec", TritiumFSpectrum, fMin / 1e9, fMax / 1e9);
  fFSpec.SetLineWidth(3);
  fFSpec.SetTitle("; f [GHz]; d#Gamma/dE [s^{-1} eV^{-1}]");
  fFSpec.SetNpx(500);

  fout.cd();
  fESpec.Write();
  fFSpec.Write();
  fout.Close();

  const double decayRateInt{fESpec.Integral(eMin, eMax)};  // s^-1 atom^-1
  const double nAtomsObs{1e20};                            // Number of atoms
  const double seconds1Yr{60 * 60 * 24 * 365};
  const double nDecays1Yr{decayRateInt * nAtomsObs * seconds1Yr};
  std::cout << "Decay rate = " << decayRateInt
            << " s^-1 atom^-1\tN decays in 1 year = " << nDecays1Yr
            << std::endl;

  // Try and do a calculation for CRESDA
  const double densCRESDA{1e18}; // m^-3
  const double volCRESDA{TMath::Pi() * 80e-3 * 30e-3 * 30e-3}; // Uniform field region
  const double nAtomsCRESDA{densCRESDA * volCRESDA};
  const double nDecaysCRESDA{nAtomsCRESDA * decayRateInt * seconds1Yr};
  std::cout << "CRESDA0 volume is " << volCRESDA << " m^3, " << volCRESDA * 1e6 << " cm^3\n";
  std::cout << "Number of atoms in CRESDA0 = " << nAtomsCRESDA << std::endl;
  std::cout << "Visible CRESDA decays in one year = " << nDecaysCRESDA << std::endl;
  return 0;
}