/*
  SignalQUick.cxx
*/

#include "SignalProcessing/SignalQuick.h"

#include <iostream>

#include "BasicFunctions/BasicFunctions.h"
#include "BasicFunctions/EMFunctions.h"
#include "TFile.h"
#include "TMath.h"
#include "TTreeReader.h"
#include "TTreeReaderValue.h"
#include "TVector3.h"

rad::SignalQuick::SignalQuick(TString trajectoryFilePath,
                              std::shared_ptr<IAntenna> ant, LocalOscillator lo,
                              double sRate,
                              std::vector<GaussianNoise> noiseTerms)
    : localOsc(lo), antenna(ant) {
  grVITime.GetYaxis()->SetTitle("V_{I}");
  grVQTime.GetYaxis()->SetTitle("V_{Q}");
  grVITime.GetXaxis()->SetTitle("Time [s]");
  grVQTime.GetXaxis()->SetTitle("Time [s]");

  antennaPos = ant->GetAntennaPosition();

  // Check if input file opens properly
  SetUpTree(trajectoryFilePath);

  // Set TTree up to be read out
  double sampleTime{0};    // Second sample time in seconds
  double sample10Time{0};  // First sample time in seconds
  unsigned int sampleNum{0};
  unsigned int sample10Num{0};
  double sampleStepSize{1 / sRate};
  double sample10StepSize{1 / (10 * sRate)};

  // Create some intermediate level graphs before final sampling
  auto grVIInter = new TGraph();
  auto grVQInter = new TGraph();

  // Loop through tree entries
  // Initially we are just doing the the higher frequency sampling
  for (unsigned int iE{0}; iE < inputTree->GetEntries(); iE++) {
    inputTree->GetEntry(iE);
    AddNewTimes(time, TVector3(xPos, yPos, zPos));

    // Do we need to sample now?
    if (time >= sample10Time) {
      // Yes we do
      // First get the retarded time to calculate the fields at
      double tr{GetRetardedTime(sample10Time)};

      double vi{CalcVoltage(tr)};
      double vq{CalcVoltage(tr)};

      // Go back to correct tree entry
      inputTree->GetEntry(iE);
      vi *= localOsc.GetInPhaseComponent(sample10Time);
      vq *= localOsc.GetQuadratureComponent(sample10Time);

      grVIInter->SetPoint(grVIInter->GetN(), sample10Time, vi);
      grVQInter->SetPoint(grVQInter->GetN(), sample10Time, vq);

      sample10Num++;
      sample10Time = double(sample10Num) * sample10StepSize;
    } else {
      continue;
    }
  }

  // Can now safely close the input file
  CloseInputFile();

  std::cout << "Filtering the signal...\n";
  auto grVIFiltered = BandPassFilter(grVIInter, 0, sRate / 2);
  auto grVQFiltered = BandPassFilter(grVQInter, 0, sRate / 2);
  delete grVIInter;
  delete grVQInter;

  std::cout << "Second sampling...\n";
  for (int i{0}; i < grVIFiltered->GetN(); i++) {
    // We need to sample every 10th point
    if (i % 10 == 0) {
      grVITime.SetPoint(i, grVIFiltered->GetPointX(i),
                        grVIFiltered->GetPointY(i));
      grVQTime.SetPoint(i, grVQFiltered->GetPointX(i),
                        grVQFiltered->GetPointY(i));
    }
  }
  delete grVIFiltered;
  delete grVQFiltered;
}

double rad::SignalQuick::CalcVoltage(double tr) {
  if (tr == -1) {
    // The voltage is from before the signal has reached the antenna
    return 0;
  } else {
    inputTree->GetEntry(0);
    const double fileStartTime{time};
    inputTree->GetEntry(inputTree->GetEntries() - 1);
    const double fileEndTime{time};
    const double pntsPerTime{double(inputTree->GetEntries()) /
                             (fileEndTime - fileStartTime)};
    // We actually have to calculate the voltage
    // Start off with a first guess
    int firstGuessTInd{int(round(pntsPerTime * (tr - fileStartTime)))};

    // Find the appropriate time
    inputTree->GetEntry(firstGuessTInd);
    double firstGuessTime{time};
    unsigned int correctIndex{0};
    if (firstGuessTime == tr) {
      // Easy, no need for interpolation
      TVector3 pos(xPos, yPos, zPos);
      TVector3 vel(xVel, yVel, zVel);
      TVector3 acc(xAcc, yAcc, zAcc);
      ROOT::Math::XYZVector eField{CalcEField(antennaPos, pos, vel, acc)};
      TVector3 eField2(eField.X(), eField.Y(), eField.Z());
      double voltage{(eField2.Dot(antenna->GetETheta(pos)) +
                      eField2.Dot(antenna->GetEPhi(pos))) *
                     antenna->GetHEff()};
      voltage /= 2.0;
      return voltage;
    } else if (firstGuessTime < tr) {
      // We are searching upwards
      for (int i{firstGuessTInd}; i < inputTree->GetEntries() - 2; i++) {
        inputTree->GetEntry(i);
        double lowerPoint{time};
        inputTree->GetEntry(i + 1);
        double upperPoint{time};
        if (tr > lowerPoint && tr < upperPoint) {
          correctIndex = i;
          break;
        }
      }
    } else if (firstGuessTime > tr) {
      // We are searching downwards
      for (int i{firstGuessTInd}; i >= 0; i--) {
        inputTree->GetEntry(i);
        double lowerPoint{time};
        inputTree->GetEntry(i + 1);
        double upperPoint{time};
        if (tr > lowerPoint && tr < upperPoint) {
          correctIndex = i;
          break;
        }
      }
    }

    // We have the relevant index so we can now do some interpolation
    std::vector<double> timeVals(4);
    std::vector<double> vVals(4);
    if (correctIndex == 0) {
      timeVals.at(0) = 0;
      vVals.at(0) = 0;
    } else {
      inputTree->GetEntry(correctIndex - 1);
      timeVals.at(0) = time;
      TVector3 pos(xPos, yPos, zPos);
      TVector3 vel(xVel, yVel, zVel);
      TVector3 acc(xAcc, yAcc, zAcc);
      ROOT::Math::XYZVector eField{CalcEField(antennaPos, pos, vel, acc)};
      TVector3 eField2(eField.X(), eField.Y(), eField.Z());
      double voltage{(eField2.Dot(antenna->GetETheta(pos)) +
                      eField2.Dot(antenna->GetEPhi(pos))) *
                     antenna->GetHEff()};
      voltage /= 2.0;
      vVals.at(0) = voltage;
    }

    // Add the rest of the elements of the vectors
    for (unsigned int iEl{1}; iEl <= 3; iEl++) {
      inputTree->GetEntry(correctIndex + iEl - 1);
      timeVals.at(iEl) = time;
      TVector3 pos(xPos, yPos, zPos);
      TVector3 vel(xVel, yVel, zVel);
      TVector3 acc(xAcc, yAcc, zAcc);
      ROOT::Math::XYZVector eField{CalcEField(antennaPos, pos, vel, acc)};
      TVector3 eField2(eField.X(), eField.Y(), eField.Z());
      double voltage{(eField2.Dot(antenna->GetETheta(pos)) +
                      eField2.Dot(antenna->GetEPhi(pos))) *
                     antenna->GetHEff()};
      voltage /= 2.0;
      vVals.at(iEl) = voltage;
    }

    // Now actually do the cubic interpolation
    double vInterp{CubicInterpolation(timeVals, vVals, tr)};
    return vInterp;
  }
}

void rad::SignalQuick::AddNewTimes(double time, TVector3 ePos) {
  timeVec.push_back(time);
  // Now calculate advanced time for this point
  double ta{time + (ePos - antennaPos).Mag() / TMath::C()};
  advancedTimeVec.push_back(ta);

  // If the vectors are getting too long, get rid of the first element
  if (timeVec.size() > 10000) {
    timeVec.erase(timeVec.begin());
    advancedTimeVec.erase(advancedTimeVec.begin());
  }
}

double rad::SignalQuick::GetRetardedTime(double ts) {
  unsigned int chosenInd{0};

  // Check if the sample time is before the signal has reached the antenna
  if (ts < advancedTimeVec.at(0)) {
    return -1;
  } else {
    // Find the values which the sample time lives in between
    // Firstly get a good first guess of where to start looking
    unsigned int firstGuess{GetFirstGuessPoint(ts)};
    // Check which direction to search in
    if (advancedTimeVec.at(firstGuess) < ts) {
      // We are searching upwards
      for (size_t i{firstGuess}; i < advancedTimeVec.size() - 2; i++) {
        // Aim to find the advanced time values the ssample point lies between
        if (ts > advancedTimeVec.at(i) && ts < advancedTimeVec.at(i + 1)) {
          chosenInd = i;
          break;
        }
      }
    } else {
      // We are searching downwards
      for (size_t i{firstGuess}; i >= 0; i--) {
        if (ts > advancedTimeVec.at(i) && ts < advancedTimeVec.at(i + 1)) {
          chosenInd = i;
          break;
        }
      }
    }

    // Now we have the relevant index, we need to interpolate
    // This should be the retarded time
    std::vector<double> timeVals(4);
    std::vector<double> advancedTimeVals(4);
    if (chosenInd == 0) {
      timeVals.at(0) = 0;
      advancedTimeVals.at(0) = 0;
    } else {
      timeVals.at(0) = timeVec.at(chosenInd - 1);
      advancedTimeVals.at(0) = advancedTimeVec.at(chosenInd - 1);
    }
    timeVals.at(1) = timeVec.at(chosenInd);
    timeVals.at(2) = timeVec.at(chosenInd + 1);
    timeVals.at(3) = timeVec.at(chosenInd + 2);
    advancedTimeVals.at(1) = advancedTimeVec.at(chosenInd);
    advancedTimeVals.at(2) = advancedTimeVec.at(chosenInd + 1);
    advancedTimeVals.at(3) = advancedTimeVec.at(chosenInd + 2);
    return CubicInterpolation(advancedTimeVals, timeVals, ts);
  }
}

void rad::SignalQuick::SetUpTree(TString filePath) {
  // Open the input file
  OpenInputFile(filePath);

  inputTree = (TTree *)inputFile->Get("tree");
  inputTree->SetBranchAddress("time", &time);
  inputTree->SetBranchAddress("xPos", &xPos);
  inputTree->SetBranchAddress("yPos", &yPos);
  inputTree->SetBranchAddress("zPos", &zPos);
  inputTree->SetBranchAddress("xVel", &xVel);
  inputTree->SetBranchAddress("yVel", &yVel);
  inputTree->SetBranchAddress("zVel", &zVel);
  inputTree->SetBranchAddress("xAcc", &xAcc);
  inputTree->SetBranchAddress("yAcc", &yAcc);
  inputTree->SetBranchAddress("zAcc", &zAcc);
}

void rad::SignalQuick::OpenInputFile(TString filePath) {
  // Check if input file opens properly
  inputFile = std::make_unique<TFile>(filePath, "READ");
  if (!inputFile->IsOpen()) {
    std::cout << "Couldn't open file! Exiting.\n";
    exit(1);
  }
}

unsigned int rad::SignalQuick::GetFirstGuessPoint(double ts) {
  const size_t taVecSize{advancedTimeVec.size()};
  const double pntsPerTime{
      double(taVecSize) /
      (advancedTimeVec.at(taVecSize - 1) - advancedTimeVec.at(0))};

  size_t firstGuessPnt{
      size_t(round(pntsPerTime * (ts - advancedTimeVec.at(0))))};
  if (firstGuessPnt < 0) firstGuessPnt = 0;
  return firstGuessPnt;
}

void rad::SignalQuick::CloseInputFile() {
  delete inputTree;
  inputFile->Close();
}