// FieldClasses.cxx

#include "FieldClasses/FieldClasses.h"

#include <cassert>
#include <cmath>
#include <vector>

#include "Antennas/IAntenna.h"
#include "BasicFunctions/BasicFunctions.h"
#include "BasicFunctions/Constants.h"
#include "Math/Point3D.h"
#include "Math/Vector3D.h"
#include "TAxis.h"
#include "TFile.h"
#include "TGraph.h"
#include "TSpline.h"
#include "TTree.h"
#include "TVector3.h"

rad::FieldPoint::~FieldPoint() {
  delete EField[0];
  delete EField[1];
  delete EField[2];
  delete BField[0];
  delete BField[1];
  delete BField[2];

  delete pos[0];
  delete pos[1];
  delete pos[2];

  delete tPrime;
}

// Parametrised constructor
rad::FieldPoint::FieldPoint(TString trajectoryFilePath, IAntenna* myAnt) {
  EField[0] = new TGraph();
  EField[1] = new TGraph();
  EField[2] = new TGraph();
  BField[0] = new TGraph();
  BField[1] = new TGraph();
  BField[2] = new TGraph();
  pos[0] = new TGraph();
  pos[1] = new TGraph();
  pos[2] = new TGraph();
  tPrime = new TGraph();

  // Now check that the input file exists
  TFile* f = new TFile(trajectoryFilePath, "read");
  assert(f);
  TTree* tree = (TTree*)f->Get("tree");
  assert(tree);
  double timeTemp;
  tree->SetBranchAddress("time", &timeTemp);
  tree->GetEntry(0);
  fileStartTime = timeTemp;
  delete tree;
  f->Close();
  delete f;
  inputFile = trajectoryFilePath;

  myAntenna = myAnt;

  minCutTime = 0.0;
  maxCutTime = 0.0;
}

// Copy constructor
rad::FieldPoint::FieldPoint(const FieldPoint& fp) {
  inputFile = fp.inputFile;

  // Clone the field graphs
  for (int coord = 0; coord < 3; coord++) {
    EField[coord] = (TGraph*)fp.EField[coord]->Clone();
    BField[coord] = (TGraph*)fp.BField[coord]->Clone();
    pos[coord] = (TGraph*)fp.pos[coord]->Clone();
  }
  tPrime = (TGraph*)fp.tPrime->Clone();
  fileStartTime = fp.fileStartTime;
  myAntenna = fp.myAntenna;

  minCutTime = fp.minCutTime;
  maxCutTime = fp.maxCutTime;
}

void rad::FieldPoint::ResetFields() {
  EField[0]->Clear();
  EField[1]->Clear();
  EField[2]->Clear();
  BField[0]->Clear();
  BField[1]->Clear();
  BField[2]->Clear();
  pos[0]->Clear();
  pos[1]->Clear();
  pos[2]->Clear();
  tPrime->Clear();

  for (int coord = 0; coord < 3; coord++) {
    delete EField[coord];
    delete BField[coord];
    delete pos[coord];
  }
  delete tPrime;

  for (int coord = 0; coord < 3; coord++) {
    EField[coord] = new TGraph();
    BField[coord] = new TGraph();
    pos[coord] = new TGraph();
    tPrime = new TGraph();
  }
}

/*
  Takes the inputted TGraph and produces the same variable plotted using the
  retarded time grOriginal is the inputted TGraph returns a TGraph which takes
  the retarded time into account
 */
TGraph* rad::FieldPoint::MakeRetardedTimeGraph(const TGraph* grOriginal) {
  TSpline3* sptPrime = new TSpline3("sptPrime", tPrime);
  TSpline3* spgrOriginal = new TSpline3("spgrOriginal", grOriginal);

  // Need to work out the first time that we should take this spline from
  double timeToStart = tPrime->GetPointX(0);

  TGraph* grOut = new TGraph();
  for (int i = 0; i < grOriginal->GetN(); i++) {
    if (grOriginal->GetPointX(i) < timeToStart) continue;
    double tRet = sptPrime->Eval(grOriginal->GetPointX(i));
    grOut->SetPoint(grOut->GetN(), grOriginal->GetPointX(i),
                    spgrOriginal->Eval(tRet));
  }
  setGraphAttr(grOut);
  grOut->GetXaxis()->SetTitle("Time [s]");
  grOut->GetYaxis()->SetTitle(grOriginal->GetYaxis()->GetTitle());

  delete sptPrime;
  delete spgrOriginal;
  return grOut;
}

TGraph* rad::FieldPoint::TrimGraphToTime(const TGraph* grIn) {
  TGraph* grOut = new TGraph();
  setGraphAttr(grOut);
  grOut->GetXaxis()->SetTitle("Time [s]");
  grOut->GetYaxis()->SetTitle(grIn->GetYaxis()->GetTitle());
  for (int i = 0; i < grIn->GetN(); i++) {
    double thisTime = grIn->GetPointX(i);
    if (thisTime < minCutTime) continue;
    if (thisTime >= maxCutTime) break;

    grOut->SetPoint(grOut->GetN(), thisTime, grIn->GetPointY(i));
  }
  return grOut;
}

// From an input TFile generate the E and B fields for a given time
// maxTime is the final time in seconds (if less than the time in the file)
void rad::FieldPoint::GenerateFields(const double minTime,
                                     const double maxTime) {
  ResetFields();

  minCutTime = minTime;
  maxCutTime = maxTime;

  TFile* fin = new TFile(inputFile, "READ");
  assert(fin);
  TTree* tree = (TTree*)fin->Get("tree");

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

  ROOT::Math::XYZPoint antennaPoint((myAntenna->GetAntennaPosition()).X(),
                                    (myAntenna->GetAntennaPosition()).Y(),
                                    (myAntenna->GetAntennaPosition()).Z());

  tree->GetEntry(0);
  const double t0 = time;
  tree->GetEntry(1);
  const double t1 = time;
  const double timeStepSize = t1 - t0;

  double minGenTime,
      maxGenTime;  // Minimum and maximum time to generate the fields between
  // If we are at the start of the file then work as normal
  if (fileStartTime == minTime) {
    minGenTime = minTime;
    maxGenTime = maxTime;
  } else {
    // Generate fields a small amount of time earlier than asked for
    // This allows for the generation of retarded time graphs which link up
    // across time chunks
    minGenTime = minTime - 4e-9;
    maxGenTime = maxTime;
  }

  // Loop through the entries and get the fields at each point
  for (int e = 0; e < tree->GetEntries(); e++) {
    tree->GetEntry(e);
    if (time < minGenTime) continue;
    if (time > maxGenTime) break;

    if (std::fmod(time, 1e-6) < timeStepSize) {
      std::cout << time << " seconds generated..." << std::endl;
    }

    ROOT::Math::XYZPoint ePos(xPos, yPos, zPos);
    ROOT::Math::XYZVector eVel(xVel, yVel, zVel);
    ROOT::Math::XYZVector eAcc(xAcc, yAcc, zAcc);
    ROOT::Math::XYZVector EFieldCalc =
        CalcEField(antennaPoint, ePos, eVel, eAcc);
    ROOT::Math::XYZVector BFieldCalc =
        CalcBField(antennaPoint, ePos, eVel, eAcc);

    EField[0]->SetPoint(EField[0]->GetN(), time, EFieldCalc.X());
    EField[1]->SetPoint(EField[1]->GetN(), time, EFieldCalc.Y());
    EField[2]->SetPoint(EField[2]->GetN(), time, EFieldCalc.Z());
    BField[0]->SetPoint(BField[0]->GetN(), time, BFieldCalc.X());
    BField[1]->SetPoint(BField[1]->GetN(), time, BFieldCalc.Y());
    BField[2]->SetPoint(BField[2]->GetN(), time, BFieldCalc.Z());

    pos[0]->SetPoint(pos[0]->GetN(), time, xPos);
    pos[1]->SetPoint(pos[1]->GetN(), time, yPos);
    pos[2]->SetPoint(pos[2]->GetN(), time, zPos);

    tPrime->SetPoint(tPrime->GetN(),
                     CalcTimeFromRetardedTime(antennaPoint, ePos, time), time);
  }

  delete tree;
  fin->Close();
  delete fin;
}

TGraph* rad::FieldPoint::GetEFieldTimeDomain(Coord_t coord,
                                             const bool kUseRetardedTime,
                                             int firstPoint, int lastPoint) {
  TGraph* gr = 0;
  TGraph* grOut = new TGraph();
  if (coord == kX) {
    gr = (TGraph*)EField[0]->Clone("grEx");
    grOut->GetYaxis()->SetTitle("E_{x} [V m^{-1}]");
  } else if (coord == kY) {
    gr = (TGraph*)EField[1]->Clone("grEy");
    grOut->GetYaxis()->SetTitle("E_{y} [V m^{-1}]");
  } else if (coord == kZ) {
    gr = (TGraph*)EField[2]->Clone("grEz");
    grOut->GetYaxis()->SetTitle("E_{z} [V m^{-1}]");
  }

  setGraphAttr(grOut);
  grOut->GetXaxis()->SetTitle("Time [s]");
  if (!kUseRetardedTime) {
    if (firstPoint < 0) firstPoint = 0;
    if (lastPoint < 0) lastPoint = gr->GetN() - 1;

    for (int i = firstPoint; i <= lastPoint; i++) {
      grOut->SetPoint(grOut->GetN(), gr->GetPointX(i), gr->GetPointY(i));
    }
  } else {
    TGraph* grRet = MakeRetardedTimeGraph(gr);
    if (firstPoint < 0) firstPoint = 0;
    if (lastPoint < 0) lastPoint = grRet->GetN() - 1;

    for (int i = firstPoint; i <= lastPoint; i++) {
      grOut->SetPoint(grOut->GetN(), grRet->GetPointX(i), grRet->GetPointY(i));
    }
    delete grRet;
  }
  delete gr;
  TGraph* grTrimmed = TrimGraphToTime(grOut);
  delete grOut;
  return grTrimmed;
}

TGraph* rad::FieldPoint::GetPositionTimeDomain(Coord_t coord,
                                               const bool kUseRetardedTime,
                                               int firstPoint, int lastPoint) {
  TGraph* gr = 0;
  TGraph* grOut = new TGraph();
  if (coord == kX) {
    gr = (TGraph*)pos[0]->Clone("grPosx");
    grOut->GetYaxis()->SetTitle("x [m]");
  } else if (coord == kY) {
    gr = (TGraph*)pos[1]->Clone("grPosy");
    grOut->GetYaxis()->SetTitle("y [m]");
  } else if (coord == kZ) {
    gr = (TGraph*)pos[2]->Clone("grPosz");
    grOut->GetYaxis()->SetTitle("z [m]");
  }

  setGraphAttr(grOut);
  grOut->GetXaxis()->SetTitle("Time [s]");
  if (!kUseRetardedTime) {
    if (firstPoint < 0) firstPoint = 0;
    if (lastPoint < 0) lastPoint = gr->GetN() - 1;

    for (int i = firstPoint; i <= lastPoint; i++) {
      grOut->SetPoint(grOut->GetN(), gr->GetPointX(i), gr->GetPointY(i));
    }
  } else {
    TGraph* grRet = MakeRetardedTimeGraph(gr);
    if (firstPoint < 0) firstPoint = 0;
    if (lastPoint < 0) lastPoint = grRet->GetN() - 1;

    for (int i = firstPoint; i <= lastPoint; i++) {
      grOut->SetPoint(grOut->GetN(), grRet->GetPointX(i), grRet->GetPointY(i));
    }
    delete grRet;
  }
  delete gr;
  TGraph* grTrimmed = TrimGraphToTime(grOut);
  delete grOut;
  return grTrimmed;
}

TGraph* rad::FieldPoint::GetEFieldMagTimeDomain(const bool kUseRetardedTime) {
  TGraph* grMag = new TGraph();
  assert((EField[0]->GetN() == EField[1]->GetN()) &&
         (EField[1]->GetN() == EField[2]->GetN()));

  for (int i = 0; i < EField[0]->GetN(); i++) {
    double mag =
        sqrt(pow(EField[0]->GetPointY(i), 2) + pow(EField[1]->GetPointY(i), 2) +
             pow(EField[2]->GetPointY(i), 2));
    grMag->SetPoint(grMag->GetN(), EField[0]->GetPointX(i), mag);
  }
  setGraphAttr(grMag);
  grMag->GetXaxis()->SetTitle("Time [s]");
  grMag->GetYaxis()->SetTitle("|E| [V m^{-1}]");

  if (!kUseRetardedTime) {
    TGraph* grTrimmed = TrimGraphToTime(grMag);
    delete grMag;
    return grTrimmed;
  } else {
    TGraph* grMagRet = MakeRetardedTimeGraph(grMag);
    delete grMag;
    TGraph* grTrimmed = TrimGraphToTime(grMagRet);
    delete grMagRet;
    return grTrimmed;
  }
}

TGraph* rad::FieldPoint::GetBFieldTimeDomain(Coord_t coord,
                                             const bool kUseRetardedTime) {
  TGraph* gr = 0;
  if (coord == kX) {
    gr = (TGraph*)BField[0]->Clone("grBx");
    gr->GetYaxis()->SetTitle("B_{x} [T]");
  } else if (coord == kY) {
    gr = (TGraph*)BField[1]->Clone("grBy");
    gr->GetYaxis()->SetTitle("B_{y} [T]");
  } else if (coord == kZ) {
    gr = (TGraph*)BField[2]->Clone("grBz");
    gr->GetYaxis()->SetTitle("B_{z} [T]");
  }
  setGraphAttr(gr);
  gr->GetXaxis()->SetTitle("Time [s]");

  if (!kUseRetardedTime) {
    TGraph* grTrimmed = TrimGraphToTime(gr);
    delete gr;
    return grTrimmed;
  } else {
    TGraph* grRet = MakeRetardedTimeGraph(gr);
    delete gr;
    TGraph* grTrimmed = TrimGraphToTime(grRet);
    delete grRet;
    return grTrimmed;
  }
}

TGraph* rad::FieldPoint::GetBFieldMagTimeDomain(const bool kUseRetardedTime) {
  TGraph* grMag = new TGraph();
  assert((BField[0]->GetN() == BField[1]->GetN()) &&
         (BField[1]->GetN() == BField[2]->GetN()));

  for (int i = 0; i < BField[0]->GetN(); i++) {
    double mag =
        sqrt(pow(BField[0]->GetPointY(i), 2) + pow(BField[1]->GetPointY(i), 2) +
             pow(BField[2]->GetPointY(i), 2));
    grMag->SetPoint(grMag->GetN(), BField[0]->GetPointX(i), mag);
  }
  setGraphAttr(grMag);
  grMag->GetXaxis()->SetTitle("Time [s]");
  grMag->GetYaxis()->SetTitle("|B| [T]");

  if (!kUseRetardedTime) {
    TGraph* grTrimmed = TrimGraphToTime(grMag);
    delete grMag;
    return grTrimmed;
  } else {
    TGraph* grMagRet = MakeRetardedTimeGraph(grMag);
    delete grMag;
    TGraph* grTrimmed = TrimGraphToTime(grMagRet);
    delete grMagRet;
    return grTrimmed;
  }
}

TGraph* rad::FieldPoint::GetPoyntingVecTimeDomain(Coord_t coord,
                                                  const bool kUseRetardedTime) {
  TGraph* grS = new TGraph();
  for (int i = 0; i < EField[0]->GetN(); i++) {
    double comp = 0;
    if (coord == kX) {
      grS->GetYaxis()->SetTitle("S_{x} [W m^{-2}]");
      comp = EField[1]->GetPointY(i) * BField[2]->GetPointY(i) -
             EField[2]->GetPointY(i) * BField[1]->GetPointY(i);
    } else if (coord == kY) {
      comp = EField[2]->GetPointY(i) * BField[0]->GetPointY(i) -
             EField[0]->GetPointY(i) * BField[2]->GetPointY(i);
      grS->GetYaxis()->SetTitle("S_{y} [W m^{-2}]");
    } else if (coord == kZ) {
      comp = EField[0]->GetPointY(i) * BField[1]->GetPointY(i) -
             EField[1]->GetPointY(i) * BField[0]->GetPointY(i);
      grS->GetYaxis()->SetTitle("S_{z} [W m^{-2}]");
    }
    comp /= MU0;
    grS->SetPoint(i, EField[0]->GetPointX(i), comp);
  }
  setGraphAttr(grS);
  grS->GetXaxis()->SetTitle("Time [s]");

  if (!kUseRetardedTime) {
    TGraph* grTrimmed = TrimGraphToTime(grS);
    delete grS;
    return grTrimmed;
  } else {
    TGraph* grSRet = MakeRetardedTimeGraph(grS);
    delete grS;
    TGraph* grTrimmed = TrimGraphToTime(grSRet);
    delete grSRet;
    return grTrimmed;
  }
}

TGraph* rad::FieldPoint::GetPoyntingMagTimeDomain(const bool kUseRetardedTime) {
  TGraph* grSMag = new TGraph();
  TGraph* grEMag = GetEFieldMagTimeDomain();
  TGraph* grBMag = GetBFieldMagTimeDomain();
  assert(grEMag->GetN() == grBMag->GetN());

  for (int i = 0; i < grEMag->GetN(); i++) {
    double smag = grEMag->GetPointY(i) * grBMag->GetPointY(i) / MU0;
    grSMag->SetPoint(grSMag->GetN(), grEMag->GetPointX(i), smag);
  }
  delete grEMag;
  delete grBMag;

  setGraphAttr(grSMag);
  grSMag->GetXaxis()->SetTitle("Time [s]");
  grSMag->GetYaxis()->SetTitle("|S| [W m^{-2}]");

  if (!kUseRetardedTime) {
    TGraph* grTrimmed = TrimGraphToTime(grSMag);
    delete grSMag;
    return grTrimmed;
  } else {
    TGraph* grSMagRet = MakeRetardedTimeGraph(grSMag);
    delete grSMag;
    TGraph* grTrimmed = TrimGraphToTime(grSMagRet);
    delete grSMagRet;
    return grTrimmed;
  }
}

TGraph* rad::FieldPoint::GetAntennaLoadVoltageTimeDomain(
    const bool kUseRetardedTime, int firstPoint, int lastPoint) {
  TGraph* grEx =
      GetEFieldTimeDomain(kX, kUseRetardedTime, firstPoint, lastPoint);
  TGraph* grEy =
      GetEFieldTimeDomain(kY, kUseRetardedTime, firstPoint, lastPoint);
  TGraph* grEz =
      GetEFieldTimeDomain(kZ, kUseRetardedTime, firstPoint, lastPoint);

  TGraph* grPosx =
      GetPositionTimeDomain(kX, kUseRetardedTime, firstPoint, lastPoint);
  TGraph* grPosy =
      GetPositionTimeDomain(kY, kUseRetardedTime, firstPoint, lastPoint);
  TGraph* grPosz =
      GetPositionTimeDomain(kZ, kUseRetardedTime, firstPoint, lastPoint);

  TGraph* gr = new TGraph();
  gr->GetXaxis()->SetTitle("Time [s]");
  gr->GetYaxis()->SetTitle("Voltage [V]");
  setGraphAttr(gr);

  for (int i = 0; i < grEx->GetN(); i++) {
    TVector3 EField(grEx->GetPointY(i), grEy->GetPointY(i), grEz->GetPointY(i));
    TVector3 ePos(grPosx->GetPointY(i), grPosy->GetPointY(i),
                  grPosz->GetPointY(i));
    double voltage = (EField.Dot(myAntenna->GetETheta(ePos)) +
                      EField.Dot(myAntenna->GetEPhi(ePos))) *
                     myAntenna->GetHEff();
    voltage /= 2.0;  // Account for re-radiated power
    gr->SetPoint(gr->GetN(), grEx->GetPointX(i), voltage);
  }

  delete grEx;
  delete grEy;
  delete grEz;
  delete grPosx;
  delete grPosy;
  delete grPosz;

  return gr;
}

TGraph* rad::FieldPoint::GetAntennaPowerTimeDomain(bool kUseRetardedTime) {
  // First of all get the Poynting vector at the point
  TGraph* grS{GetPoyntingMagTimeDomain(kUseRetardedTime)};
  // Get the graphs of the electron position
  TGraph* grX{GetPositionTimeDomain(kX, kUseRetardedTime)};
  TGraph* grY{GetPositionTimeDomain(kY, kUseRetardedTime)};
  TGraph* grZ{GetPositionTimeDomain(kZ, kUseRetardedTime)};

  TGraph* grPower = new TGraph();
  setGraphAttr(grPower);
  grPower->GetYaxis()->SetTitle("Collected power [W]");
  grPower->GetXaxis()->SetTitle("Time [s]");
  for (int n{0}; n < grS->GetN(); n++) {
    TVector3 ePos(grX->GetPointY(n), grY->GetPointY(n), grZ->GetPointY(n));
    double AEff{myAntenna->GetAEff(ePos)};
    grPower->SetPoint(n, grS->GetPointX(n), grS->GetPointY(n) * AEff);
  }

  delete grS;
  delete grX;
  delete grY;
  delete grZ;
  return grPower;
}

TGraph* rad::FieldPoint::GetAntennaLoadPowerTimeDomain(
    const double loadResistance, const bool kUseRetardedTime, int firstPoint,
    int lastPoint) {
  TGraph* gr =
      GetAntennaLoadVoltageTimeDomain(kUseRetardedTime, firstPoint, lastPoint);
  for (int i = 0; i < gr->GetN(); i++) {
    gr->SetPointY(i, gr->GetPointY(i) * gr->GetPointY(i) / loadResistance);
  }
  return gr;
}

////////////////////////////////////////////////////////////
/////////////// Frequency domain functions /////////////////
////////////////////////////////////////////////////////////

TGraph* rad::FieldPoint::GetEFieldPeriodogram(Coord_t coord,
                                              const bool kUseRetardedTime) {
  TGraph* grIn = GetEFieldTimeDomain(coord, kUseRetardedTime);
  TGraph* grFFT = MakePowerSpectrumPeriodogram(grIn);

  if (coord == kX) {
    grFFT->GetYaxis()->SetTitle("E_{x}^{2} [V^{2} m^{-2}]");
  } else if (coord == kY) {
    grFFT->GetYaxis()->SetTitle("E_{y}^{2} [V^{2} m^{-2}]");
  } else if (coord == kZ) {
    grFFT->GetYaxis()->SetTitle("E_{z}^{2} [V^{2} m^{-2}]");
  }
  setGraphAttr(grFFT);
  grFFT->GetXaxis()->SetTitle("Frequency [Hz]");

  delete grIn;
  return grFFT;
}

TGraph* rad::FieldPoint::GetEFieldPowerSpectrumNorm(
    Coord_t coord, const bool kUseRetardedTime) {
  TGraph* grIn = GetEFieldTimeDomain(coord, kUseRetardedTime);
  TGraph* grFFT = MakePowerSpectrumNorm(grIn);

  if (coord == kX) {
    grFFT->GetYaxis()->SetTitle("E_{x}^{2} (#Deltat)^{2} [V^{2} m^{-2} s^{2}]");
  } else if (coord == kY) {
    grFFT->GetYaxis()->SetTitle("E_{y}^{2} (#Deltat)^{2} [V^{2} m^{-2} s^{2}]");
  } else if (coord == kZ) {
    grFFT->GetYaxis()->SetTitle("E_{z}^{2} (#Deltat)^{2} [V^{2} m^{-2} s^{2}]");
  }
  setGraphAttr(grFFT);
  grFFT->GetXaxis()->SetTitle("Frequency [Hz]");

  delete grIn;
  return grFFT;
}

// Total electric field in the frequency domain
TGraph* rad::FieldPoint::GetTotalEFieldPeriodogram(
    const bool kUseRetardedTime) {
  TGraph* grTotal = new TGraph();
  TGraph* grX = GetEFieldPeriodogram(kX, kUseRetardedTime);
  TGraph* grY = GetEFieldPeriodogram(kY, kUseRetardedTime);
  TGraph* grZ = GetEFieldPeriodogram(kZ, kUseRetardedTime);

  for (int n = 0; n < grX->GetN(); n++) {
    double tot = grX->GetPointY(n) + grY->GetPointY(n) + grZ->GetPointY(n);
    grTotal->SetPoint(n, grX->GetPointX(n), tot);
  }
  setGraphAttr(grTotal);
  grTotal->GetXaxis()->SetTitle("Frequency [Hz]");

  return grTotal;
}

TGraph* rad::FieldPoint::GetTotalEFieldPowerSpectrumNorm(
    const bool kUseRetardedTime) {
  TGraph* grTotal = new TGraph();
  TGraph* grX = GetEFieldPowerSpectrumNorm(kX, kUseRetardedTime);
  TGraph* grY = GetEFieldPowerSpectrumNorm(kY, kUseRetardedTime);
  TGraph* grZ = GetEFieldPowerSpectrumNorm(kZ, kUseRetardedTime);

  for (int n = 0; n < grX->GetN(); n++) {
    double tot = grX->GetPointY(n) + grY->GetPointY(n) + grZ->GetPointY(n);
    grTotal->SetPoint(n, grX->GetPointX(n), tot);
  }
  setGraphAttr(grTotal);
  grTotal->GetYaxis()->SetTitle("E^{2} (#Deltat)^{2} [V^{2} m^{-2} s^{2}]");
  grTotal->GetXaxis()->SetTitle("Frequency [Hz]");

  delete grX;
  delete grY;
  delete grZ;
  return grTotal;
}

TGraph* rad::FieldPoint::GetAntennaLoadPowerSpectrumNorm(
    const double resistance, const bool kUseRetardedTime, int firstPoint,
    int lastPoint) {
  TGraph* grVoltage =
      GetAntennaLoadVoltageTimeDomain(kUseRetardedTime, firstPoint, lastPoint);
  TGraph* grPower = MakePowerSpectrumNorm(grVoltage);

  for (int i = 0; i < grPower->GetN(); i++) {
    double powerWatts = grPower->GetPointY(i) / resistance;
    grPower->SetPoint(i, grPower->GetPointX(i), powerWatts);
  }
  setGraphAttr(grPower);
  grPower->GetXaxis()->SetTitle("Frequency [Hz]");
  grPower->GetYaxis()->SetTitle("Power #times (#Delta t)^{2} [W s^{2}]");

  delete grVoltage;
  return grPower;
}

// Assorted useful functions
double rad::FieldPoint::GetFinalTime() {
  TFile* fin = new TFile(inputFile, "READ");
  assert(fin);
  TTree* tree = (TTree*)fin->Get("tree");
  double lastTime;
  tree->SetBranchAddress("time", &lastTime);
  tree->GetEntry(tree->GetEntries() - 1);
  delete tree;
  fin->Close();
  delete fin;
  return lastTime;
}

double rad::FieldPoint::GetSampleRate() {
  TFile* fin = new TFile(inputFile, "READ");
  assert(fin);
  TTree* tree = (TTree*)fin->Get("tree");
  double time;
  tree->SetBranchAddress("time", &time);
  double time0, time1;
  tree->GetEntry(0);
  time0 = time;
  tree->GetEntry(1);
  time1 = time;
  double fs = 1.0 / (time1 - time0);
  delete tree;
  fin->Close();
  delete fin;
  return fs;
}
