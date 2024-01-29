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
#include "TFile.h"
#include "TString.h"
#include "TTree.h"
#include "TVector3.h"

using namespace rad;

TVector3 EFieldMagneticDipole(double rho, double theta, double t, double f,
                              double phase = 0, double m0 = 1) {
  const double omega{2 * M_PI * f};
  const double prefac{MU0 * m0 * sin(theta) / (4 * M_PI * rho)};
  const double ePhi{(omega * omega / TMath::C()) *
                        cos(omega * (t - rho / TMath::C())) +
                    (omega / rho) * sin(omega * (t - rho / TMath::C()))};
  return prefac * TVector3(-ePhi * sin(theta), ePhi * cos(theta), 0);
}

TVector3 EFieldMagneticDipole(TVector3 vec, double t, double f,
                              double phase = 0, double m0 = 1) {
  return EFieldMagneticDipole(vec.Perp(), vec.Theta(), t, f, phase, m0);
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
  for (int iE{0}; iE < tr->GetEntries(); iE++) {
    tr->GetEntry(iE);
    TVector3 E{(EFieldMagneticDipole(fieldPoint, time, fCyc, 0, m0))};
    grMagDX->SetPoint(grMagDX->GetN(), time, E.X());
    grMagDY->SetPoint(grMagDY->GetN(), time, E.Y());
    grMagDZ->SetPoint(grMagDZ->GetN(), time, E.Z());
  }
  fin.Close();

  fout.cd();
  grMagDX->Write("grMagDX");
  grMagDY->Write("grMagDY");
  grMagDZ->Write("grMagDZ");

  fout.Close();
  return 0;
}