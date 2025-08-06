// InducedVoltage.cxx

#include "physics/SignalProcessing/InducedVoltage.h"
#include "physics/Antennas/IAntenna.h"
#include "physics/FieldClasses/FieldClasses.h"
#include "utilities/BasicFunctions/BasicFunctions.h"

#include "TString.h"
#include "TGraph.h"
#include "TFile.h"
#include "TTree.h"
#include "TAxis.h"
#include "TSpline.h"

#include <iostream>

rad::InducedVoltage::~InducedVoltage() {
  delete grVoltage;
  theAntennas.clear();
}

TGraph* rad::InducedVoltage::DelayVoltage(TGraph* grIn, IAntenna* ant)
{
  double delay = ant->GetTimeDelay();
  TGraph* grOut = new TGraph();
  setGraphAttr(grOut);
  TSpline3* sp = new TSpline3("sp", grIn);
  for (int iPnt = 0; iPnt < grIn->GetN(); iPnt++) {
    if (grIn->GetPointX(iPnt) - delay < grIn->GetPointX(0)) continue;
    double theTime = grIn->GetPointX(iPnt);
    grOut->SetPoint(grOut->GetN(), theTime, sp->Eval(theTime - delay));
  }
  
  delete sp;
  return grOut;
}

rad::InducedVoltage::InducedVoltage(TString trajectoryFilePath, IAntenna* myAntenna,
				    const bool kUseRetardedTime) {
  theFile = trajectoryFilePath;
  theAntennas.push_back(myAntenna);
  UseRetardedTime = kUseRetardedTime;
  grVoltage = new TGraph();
  grVoltage->GetXaxis()->SetTitle("Time [s]");
  grVoltage->GetYaxis()->SetTitle("Voltage [V]");
  
  // Get the time spacing in the input file
  TFile* file1 = new TFile(theFile, "READ");
  assert(file1);
  TTree* tree1 = (TTree*)file1->Get("tree");
  assert(tree1);
  double time;
  tree1->SetBranchAddress("time", &time);
  tree1->GetEntry(0);
  double time0 = time;
  tree1->GetEntry(1);
  double time1 = time;
  delete tree1;
  file1->Close();
  delete file1;
  const double timeStep = time1 - time0;

  const double chunkRatio = 8333333.0; // Number of points that have been determined to work
  chunkSize = chunkRatio * timeStep; // Adaptive time chunk size
}

rad::InducedVoltage::InducedVoltage(TString trajectoryFilePath, std::vector<IAntenna*> antennaVec,
				    const bool kUseRetardedTime)
{
  theFile = trajectoryFilePath;
  theAntennas = antennaVec;
  UseRetardedTime = kUseRetardedTime;
  grVoltage = new TGraph();
  grVoltage->GetXaxis()->SetTitle("Time [s]");
  grVoltage->GetYaxis()->SetTitle("Voltage [V]");
  
  // Get the time spacing in the input file
  TFile* file1 = new TFile(theFile, "READ");
  assert(file1);
  TTree* tree1 = (TTree*)file1->Get("tree");
  assert(tree1);
  double time;
  tree1->SetBranchAddress("time", &time);
  tree1->GetEntry(0);
  double time0 = time;
  tree1->GetEntry(1);
  double time1 = time;
  delete tree1;
  file1->Close();
  delete file1;
  const double timeStep = time1 - time0;

  const double chunkRatio = 8333333.0; // Number of points that have been determined to work
  chunkSize = chunkRatio * timeStep; // Adaptive time chunk size
}

void rad::InducedVoltage::ProcessTimeChunk(FieldPoint fp, double firstTime, double lastTime, double minTime, double &latestStartTime, bool kFirstAntenna)
{
  double timeDelay = fp.GetAntenna()->GetTimeDelay();
  fp.GenerateFields(firstTime-timeDelay, lastTime);

  TGraph* voltageTemp = 0;
  if (timeDelay != 0.0) {
    TGraph* voltageUnshifted = fp.GetAntennaLoadVoltageTimeDomain(UseRetardedTime);
    voltageTemp = DelayVoltage(voltageUnshifted, fp.GetAntenna());
    delete voltageUnshifted;
  }
  else {
    voltageTemp = fp.GetAntennaLoadVoltageTimeDomain(UseRetardedTime);
  }

  // Now write this to the main voltage graph
  std::cout<<"Writing to main voltage graph"<<std::endl;
  if (kFirstAntenna) {
    // This is the first antenna so write to the main graph normally
    for (int i = 0; i < voltageTemp->GetN(); i++) {
      grVoltage->SetPoint(grVoltage->GetN(), voltageTemp->GetPointX(i), voltageTemp->GetPointY(i));
    }
    latestStartTime = grVoltage->GetPointX(0);
  }
  else {
    // We have already written to this graph once
    // Need to find where to start writing this graph to
    // Is this the first time chunk in the sequence?
    if (firstTime == minTime) {
      // This is the first time chunk
      // Check if this voltage has a later start time than current limit
      if (voltageTemp->GetPointX(0) > latestStartTime)  latestStartTime = voltageTemp->GetPointX(0);

      // Now figure out where to start adding these points to the existing graph
      int startPntTmp = -1;
      double startTimeTmp = 0.0;
      // Loop through points of temporary graph to determine the start time / point
      for (int iTemp = 0; iTemp < voltageTemp->GetN(); iTemp++) {
	if (voltageTemp->GetPointX(iTemp) < latestStartTime) {
	  continue;
	}
	else {
	  startPntTmp = iTemp;
	  startTimeTmp = voltageTemp->GetPointX(iTemp);
	  break;
	}	  
      }

      // Now find the corresponding point on the main graph that corresponds to this time
      int startPntMain = -1;
      for (int iMain = 0; iMain < grVoltage->GetN(); iMain++) {
	if (grVoltage->GetPointX(iMain) == startTimeTmp) {
	  startPntMain = iMain;
	  break;
	}
	    
	if (iMain == grVoltage->GetN()-1) {
	  std::cout<<"We seem to have not found a matching time point. Exiting..."<<std::endl;
	  exit(1);
	}
      }

      // Now we have these two points, can add voltages to existing graph
      for (int i = 0; i < voltageTemp->GetN() - startPntTmp; i++) {
	double existingVoltage = grVoltage->GetPointY(startPntMain+i);
	grVoltage->SetPointY(startPntMain+i, existingVoltage + voltageTemp->GetPointY(i+startPntTmp));
      } 
    }
    else {
      // This is NOT the first time chunk
      // Therefore the first point of the temporary graph should match with a main graph point
      int startPntMain = -1;
      for (int iMain = 0; iMain < grVoltage->GetN(); iMain++) {
	if (grVoltage->GetPointX(iMain) == voltageTemp->GetPointX(0)) {
	  startPntMain = iMain;
	  break;
	}

	if (iMain == grVoltage->GetN()-1) {
	  std::cout<<"We seem to have not found a matching time point. Exiting..."<<std::endl;
	  exit(1);
	}
      } // Loop over main graph points

      // Now write to the main voltage graph
      for (int i = 0; i < voltageTemp->GetN(); i++) {
	double existingVoltage = grVoltage->GetPointY(startPntMain+i);
	grVoltage->SetPointY(startPntMain+i, existingVoltage + voltageTemp->GetPointY(i));
      } // Write to main voltage graph
	  
    } // This is not the first time chunk
  } // Already written to the main graph at least once
  
  delete voltageTemp;
}

void rad::InducedVoltage::GenerateVoltage(double minTime, double maxTime) {
  double latestStartTime = -DBL_MAX;
  
  // Loop over the inputted antennas
  for (int iAnt = 0; iAnt < theAntennas.size(); iAnt++) {
    FieldPoint fp(theFile, theAntennas[iAnt]);
    if (minTime == -1) minTime = 0.0;
    if (maxTime == -1) maxTime = fp.GetFinalTime();
    bool firstAntenna = (iAnt == 0);
    
    // To avoid running out of memory, generate the fields in more manageable chunks
    // Avoids having massive versions of unnecessary graphs
    double thisChunk = minTime + chunkSize;
    if (thisChunk > maxTime) thisChunk = maxTime;
    double lastChunk = minTime;

    double timeDelay = theAntennas[iAnt]->GetTimeDelay();
    std::cout<<"Generating voltages"<<std::endl;
    while (thisChunk <= maxTime && thisChunk != lastChunk) {
      ProcessTimeChunk(fp, lastChunk, thisChunk, minTime, latestStartTime, firstAntenna);
      lastChunk = thisChunk;
      thisChunk += chunkSize;
      if (thisChunk > maxTime) thisChunk = maxTime;
    } // Keep processing chunks
  } // Loop over antenna points

  std::cout<<"Removing points up to "<<latestStartTime<<std::endl;
  // Now remove the initial points that are before the first matching start time
  for (int i = 0; i < grVoltage->GetN(); i++) {
    if (grVoltage->GetPointX(0) >= latestStartTime)
      break;
    else
      grVoltage->RemovePoint(0);
  } // Remove unmatched points
}

rad::InducedVoltage::InducedVoltage(const InducedVoltage &iv) {
  grVoltage = (TGraph*)iv.grVoltage->Clone();
  theAntennas = iv.theAntennas;
  theFile = iv.theFile;
  UseRetardedTime = iv.UseRetardedTime;
  chunkSize = iv.chunkSize;
}

TGraph* rad::InducedVoltage::GetVoltageGraph() {
  TGraph* grOut = (TGraph*)grVoltage->Clone();
  setGraphAttr(grOut);
  return grOut;
}

void rad::InducedVoltage::ResetVoltage() {
  delete grVoltage;
  grVoltage = new TGraph();
}

double rad::InducedVoltage::GetFinalTime() {
  TFile *fin = new TFile(theFile, "READ");
  assert(fin);
  TTree* tree = (TTree*)fin->Get("tree");
  double lastTime;
  tree->SetBranchAddress("time", &lastTime);
  tree->GetEntry(tree->GetEntries()-1);
  delete tree;
  fin->Close();                                                                 
  delete fin;
  return lastTime;
}

double rad::InducedVoltage::GetUpperAntennaBandwidth() {
  return theAntennas[0]->GetBandwidthUpperLimit();
}

double rad::InducedVoltage::GetLowerAntennaBandwidth() {
  return theAntennas[0]->GetBandwidthLowerLimit();
}

void rad::InducedVoltage::ApplyAntennaBandwidth() {
  grVoltage = BandPassFilter(grVoltage, theAntennas[0]->GetBandwidthLowerLimit(), theAntennas[0]->GetBandwidthUpperLimit());
}

TGraph* rad::InducedVoltage::GetPowerPeriodogram(const double loadResistance) {
  std::cout<<"Creating power spectrum"<<std::endl;
  TGraph* grV = GetVoltageGraph();
  TGraph* pgram = MakePowerSpectrumPeriodogram(grV);
  ScaleGraph(pgram, 1.0 / loadResistance);
  pgram->GetXaxis()->SetTitle("Frequency [Hz]");
  pgram->GetYaxis()->SetTitle("Power [W]");
  delete grV;
  return pgram;
}
