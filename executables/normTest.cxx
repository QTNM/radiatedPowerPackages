// normTest.cxx

// STL
#include <cmath>
#include <iostream>

// ROOT includes
#include "TAxis.h"
#include "TFile.h"
#include "TGraph.h"
#include "TMath.h"
#include "TTree.h"

// My files
#include "BasicFunctions/BasicFunctions.h"
#include "FieldClasses/FieldClasses.h"

using namespace rad;

int main() {
  TFile* fin =
      new TFile("/home/sjones/work/qtnm/trajectories/90DegOnAxis.root", "read");
  TTree* tree = (TTree*)fin->Get("tree");
  double time;
  double xPos, yPos, zPos;
  double xVel, yVel, zVel;
  double xAcc, yAcc, zAcc;
  tree->SetBranchAddress("time", &time);
  tree->SetBranchAddress("xPos", &xPos);
  tree->SetBranchAddress("yPos", &yPos);
  tree->SetBranchAddress("zPos", &zPos);
  tree->SetBranchAddress("xVel", &xVel);
  tree->SetBranchAddress("yVel", &yVel);
  tree->SetBranchAddress("zVel", &zVel);
  tree->SetBranchAddress("xAcc", &xAcc);
  tree->SetBranchAddress("yAcc", &yAcc);
  tree->SetBranchAddress("zAcc", &zAcc);

  const double maxTime = 1e-6;
  TGraph* grEy = new TGraph();

  ROOT::Math::XYZPoint antennaPoint(0.02, 0, 0);

  // Loop through the entries and get the fields at each point
  for (int e = 0; e < tree->GetEntries(); e++) {
    tree->GetEntry(e);
    if (time > maxTime) break;
    ROOT::Math::XYZPoint ePos(xPos, yPos, zPos);
    ROOT::Math::XYZVector eVel(xVel, yVel, zVel);
    ROOT::Math::XYZVector eAcc(xAcc, yAcc, zAcc);
    ROOT::Math::XYZVector EFieldCalc =
        CalcEField(antennaPoint, ePos, eVel, eAcc);
    ROOT::Math::XYZVector BFieldCalc =
        CalcBField(antennaPoint, ePos, eVel, eAcc);
    grEy->SetPoint(grEy->GetN(), time, EFieldCalc.Y());
  }
  TGraph* grEyPower = MakePowerSpectrumPeriodogram(grEy);

  double vsquared = SumVoltageSquared(grEy, -1, -1) / grEy->GetN();
  double powerSum = SumPower(grEyPower);
  std::cout << "vsquared, power = " << vsquared << ", " << powerSum
            << std::endl;

  return 0;
}
