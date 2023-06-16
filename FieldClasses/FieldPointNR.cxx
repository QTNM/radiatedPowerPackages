#include "FieldClasses/FieldPointNR.h"

#include "TTree.h"

// From an input TFile generate the E and B fields for a given time
// maxTime is the final time in seconds (if less than the time in the file)
void rad::FieldPointNR::GenerateFields(const double minTime, const double maxTime) {
  ResetFields();

  minCutTime = minTime;
  maxCutTime = maxTime;
  
  TFile *fin = new TFile(inputFile, "READ");
  assert(fin);
  TTree *tree = (TTree*)fin->Get("tree");

  // Set variables
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

  TVector3 antennaPoint{myAntenna->GetAntennaPosition()};

  tree->GetEntry(0);
  const double t0 = time;
  tree->GetEntry(1);
  const double t1 = time;
  const double timeStepSize = t1 - t0;

  double minGenTime, maxGenTime; // Minimum and maximum time to generate the fields between
  // If we are at the start of the file then work as normal
  if (fileStartTime == minTime) {
    minGenTime = minTime;
    maxGenTime = maxTime;
  }
  else {
    // Generate fields a small amount of time earlier than asked for
    // This allows for the generation of retarded time graphs which link up across time chunks
    minGenTime = minTime - 4e-9;
    maxGenTime = maxTime;
  }
  
  // Loop through the entries and get the fields at each point
  for (int e = 0; e < tree->GetEntries(); e++) {
    tree->GetEntry(e);
    if (time < minGenTime) continue;
    if (time > maxGenTime) break;

    if (std::fmod(time, 1e-6) < timeStepSize) {
      std::cout<<time<<" seconds generated..."<<std::endl;
    }
    
    TVector3 ePos(xPos, yPos, zPos);
    TVector3 eVel(xVel, yVel, zVel);
    TVector3 eAcc(xAcc, yAcc, zAcc);
    TVector3 EFieldCalc = CalcEFieldNR(antennaPoint, ePos, eVel, eAcc);
    TVector3 BFieldCalc = CalcBFieldNR(antennaPoint, ePos, eVel, eAcc);

    EField[0]->SetPoint(EField[0]->GetN(), time, EFieldCalc.X());
    EField[1]->SetPoint(EField[1]->GetN(), time, EFieldCalc.Y());
    EField[2]->SetPoint(EField[2]->GetN(), time, EFieldCalc.Z());
    BField[0]->SetPoint(BField[0]->GetN(), time, BFieldCalc.X());
    BField[1]->SetPoint(BField[1]->GetN(), time, BFieldCalc.Y());
    BField[2]->SetPoint(BField[2]->GetN(), time, BFieldCalc.Z());

    pos[0]->SetPoint(pos[0]->GetN(), time, xPos);
    pos[1]->SetPoint(pos[1]->GetN(), time, yPos);
    pos[2]->SetPoint(pos[2]->GetN(), time, zPos);

    tPrime->SetPoint(tPrime->GetN(), CalcTimeFromRetardedTime(antennaPoint, ePos, time), time);
  }
  
  delete tree;
  fin->Close();
  delete fin;
}