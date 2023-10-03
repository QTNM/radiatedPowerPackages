/*
  Signal.cxx
*/

#include "SignalProcessing/Signal.h"

#include <iostream>

#include "BasicFunctions/BasicFunctions.h"
#include "BasicFunctions/ComplexVector3.h"
#include "BasicFunctions/EMFunctions.h"
#include "TMath.h"
#include "TTreeReader.h"
#include "TTreeReaderValue.h"
#include "TVector3.h"

rad::Signal::Signal(TString trajectoryFilePath, IAntenna* ant,
                    LocalOscillator lo, double sRate,
                    std::vector<GaussianNoise> noiseTerms, double tAcq)
    : localOsc(lo), sampleRate(sRate), noiseVec(noiseTerms) {
  antenna.push_back(ant);

  CreateVoltageGraphs();

  // Check if input file opens properly
  SetUpTree(trajectoryFilePath);

  // Set file info
  GetFileInfo();

  // Figure out where we're going to generate the signal up to
  // By default, just do the whole electron trajectory file
  if (tAcq < 0) tAcq = fileEndTime;

  double sampleTime{0};    // Second sample time in seconds
  double sample10Time{0};  // First sample time in seconds
  unsigned int sample10Num{0};
  double sample10StepSize{1 / (10 * sRate)};

  // Create just one deque
  advancedTimeVec.push_back(std::deque<double>());

  // Create some intermediate level graphs before final sampling
  auto grVIInter = new TGraph();
  auto grVQInter = new TGraph();
  auto grVIBigFiltered = new TGraph();
  auto grVQBigFiltered = new TGraph();

  // Loop through tree entries
  // Initially we are just doing the the higher frequency sampling
  double printTime{0};  // seconds
  const double printInterval{5e-6};
  for (unsigned int iE{0}; iE < inputTree->GetEntries(); iE++) {
    inputTree->GetEntry(iE);
    const double entryTime{time};

    // Check we are still within the acquisition time, stop otherwise
    if (entryTime > tAcq) break;

    if (entryTime >= printTime) {
      std::cout << printTime * 1e6 << " us signal processed...\n";
      printTime += printInterval;
    }

    AddNewTimes(entryTime, TVector3(xPos, yPos, zPos));

    // Do we need to sample now?
    if (entryTime >= sample10Time) {
      // Yes we do
      // First get the retarded time to calculate the fields at
      double tr{GetRetardedTime(sample10Time, 0)};

      double vi{CalcVoltage(tr, antenna[0])};
      double vq{vi};

      DownmixVoltages(vi, vq, sample10Time);

      grVIInter->SetPoint(grVIInter->GetN(), sample10Time, vi);
      grVQInter->SetPoint(grVQInter->GetN(), sample10Time, vq);

      sample10Num++;
      sample10Time = double(sample10Num) * sample10StepSize;

      // We want to filter this signal if it's long enough
      // Pick a good number to do FFTs with (2^n preferably)
      if (grVIInter->GetN() == 32768) {
        auto grVISmallFiltered = BandPassFilter(grVIInter, 0, sRate / 2);
        auto grVQSmallFiltered = BandPassFilter(grVQInter, 0, sRate / 2);
        // Now add these points to the existing graph
        for (int iF{0}; iF < grVISmallFiltered->GetN(); iF++) {
          grVIBigFiltered->SetPoint(grVIBigFiltered->GetN(),
                                    grVISmallFiltered->GetPointX(iF),
                                    grVISmallFiltered->GetPointY(iF));
          grVQBigFiltered->SetPoint(grVQBigFiltered->GetN(),
                                    grVQSmallFiltered->GetPointX(iF),
                                    grVQSmallFiltered->GetPointY(iF));
        }
        delete grVISmallFiltered;
        delete grVQSmallFiltered;
        // Clear the current graphs and start again
        delete grVIInter;
        delete grVQInter;
        grVIInter = new TGraph();
        grVQInter = new TGraph();
      }
    } else {
      continue;
    }
  }

  // Can now safely close the input file
  CloseInputFile();

  // One last filter maybe
  if (grVIInter->GetN() > 0) {
    TGraph* grVISmallFiltered = BandPassFilter(grVIInter, 0, sRate / 2);
    TGraph* grVQSmallFiltered = BandPassFilter(grVQInter, 0, sRate / 2);
    delete grVIInter;
    delete grVQInter;
    // Add these remaining points to the main filtered graph
    for (int iF{0}; iF < grVISmallFiltered->GetN(); iF++) {
      grVIBigFiltered->SetPoint(grVIBigFiltered->GetN(),
                                grVISmallFiltered->GetPointX(iF),
                                grVISmallFiltered->GetPointY(iF));
      grVQBigFiltered->SetPoint(grVQBigFiltered->GetN(),
                                grVQSmallFiltered->GetPointX(iF),
                                grVQSmallFiltered->GetPointY(iF));
    }
    delete grVISmallFiltered;
    delete grVQSmallFiltered;
  }

  std::cout << "Second sampling...\n";
  for (int i{0}; i < grVIBigFiltered->GetN(); i++) {
    // We need to sample every 10th point
    if (i % 10 == 0) {
      grVITime->SetPoint(grVITime->GetN(), grVIBigFiltered->GetPointX(i),
                         grVIBigFiltered->GetPointY(i));
      grVQTime->SetPoint(grVQTime->GetN(), grVQBigFiltered->GetPointX(i),
                         grVQBigFiltered->GetPointY(i));
    }
  }
  delete grVIBigFiltered;
  delete grVQBigFiltered;

  // Now need to add noise (if noise terms exist)
  if (!noiseTerms.empty()) {
    std::cout << "Adding noise...\n";
    AddNoise();
  }
}

rad::Signal::Signal(TString trajectoryFilePath, std::vector<IAntenna*> ant,
                    LocalOscillator lo, double sRate,
                    std::vector<GaussianNoise> noiseTerms, double tAcq)
    : localOsc(lo), sampleRate(sRate), noiseVec(noiseTerms), antenna(ant) {
  CreateVoltageGraphs();

  // Check if input file opens properly
  SetUpTree(trajectoryFilePath);

  // Set file info
  GetFileInfo();

  double sampleTime{0};    // Second sample time in seconds
  double sample10Time{0};  // First sample time in seconds
  unsigned int sample10Num{0};
  double sample10StepSize{1 / (10 * sRate)};

  // Create number of deques equal to number of antennas
  for (size_t i{0}; i < antenna.size(); i++) {
    advancedTimeVec.push_back(std::deque<double>());
  }

  // Create some intermediate level graphs before final sampling
  auto grVIInter = new TGraph();
  auto grVQInter = new TGraph();
  auto grVIBigFiltered = new TGraph();
  auto grVQBigFiltered = new TGraph();

  // Loop through tree entries
  // Initially we are just doing the the higher frequency sampling
  double printTime{0};  // seconds
  double printInterval{5e-6};
  for (unsigned int iE{0}; iE < inputTree->GetEntries(); iE++) {
    inputTree->GetEntry(iE);
    const double entryTime{time};
    if (entryTime >= printTime) {
      std::cout << printTime * 1e6 << " us signal processed...\n";
      printTime += printInterval;
    }

    AddNewTimes(entryTime, TVector3(xPos, yPos, zPos));

    // Do we need to sample now?
    if (entryTime >= sample10Time) {
      // Yes we do
      // Loop through antennas and calculate a voltage at each
      double vi{0};
      double vq{0};
      for (size_t iAnt{0}; iAnt < antenna.size(); iAnt++) {
        // Yes we do
        // First get the retarded time to calculate the fields at
        double tr{GetRetardedTime(sample10Time, iAnt)};
        double v{CalcVoltage(tr, antenna[iAnt])};
        vi += v;
        vq += v;
      }
      DownmixVoltages(vi, vq, sample10Time);

      grVIInter->SetPoint(grVIInter->GetN(), sample10Time, vi);
      grVQInter->SetPoint(grVQInter->GetN(), sample10Time, vq);

      sample10Num++;
      sample10Time = double(sample10Num) * sample10StepSize;

      // We want to filter this signal if it's long enough
      // Pick a good number to do FFTs with (2^n preferably)
      if (grVIInter->GetN() == 32768) {
        auto grVISmallFiltered = BandPassFilter(grVIInter, 0, sRate / 2);
        auto grVQSmallFiltered = BandPassFilter(grVQInter, 0, sRate / 2);
        // Now add these points to the existing graph
        for (int iF{0}; iF < grVISmallFiltered->GetN(); iF++) {
          grVIBigFiltered->SetPoint(grVIBigFiltered->GetN(),
                                    grVISmallFiltered->GetPointX(iF),
                                    grVISmallFiltered->GetPointY(iF));
          grVQBigFiltered->SetPoint(grVQBigFiltered->GetN(),
                                    grVQSmallFiltered->GetPointX(iF),
                                    grVQSmallFiltered->GetPointY(iF));
        }
        delete grVISmallFiltered;
        delete grVQSmallFiltered;
        // Clear the current graphs and start again
        delete grVIInter;
        delete grVQInter;
        grVIInter = new TGraph();
        grVQInter = new TGraph();
      }
    } else {
      continue;
    }
  }

  // Can now safely close the input file
  CloseInputFile();

  // One last filter maybe
  if (grVIInter->GetN() > 0) {
    TGraph* grVISmallFiltered = BandPassFilter(grVIInter, 0, sRate / 2);
    TGraph* grVQSmallFiltered = BandPassFilter(grVQInter, 0, sRate / 2);
    delete grVIInter;
    delete grVQInter;
    // Add these remaining points to the main filtered graph
    for (int iF{0}; iF < grVISmallFiltered->GetN(); iF++) {
      grVIBigFiltered->SetPoint(grVIBigFiltered->GetN(),
                                grVISmallFiltered->GetPointX(iF),
                                grVISmallFiltered->GetPointY(iF));
      grVQBigFiltered->SetPoint(grVQBigFiltered->GetN(),
                                grVQSmallFiltered->GetPointX(iF),
                                grVQSmallFiltered->GetPointY(iF));
    }
    delete grVISmallFiltered;
    delete grVQSmallFiltered;
  }

  std::cout << "Second sampling...\n";
  for (int i{0}; i < grVIBigFiltered->GetN(); i++) {
    // We need to sample every 10th point
    if (i % 10 == 0) {
      grVITime->SetPoint(grVITime->GetN(), grVIBigFiltered->GetPointX(i),
                         grVIBigFiltered->GetPointY(i));
      grVQTime->SetPoint(grVQTime->GetN(), grVQBigFiltered->GetPointX(i),
                         grVQBigFiltered->GetPointY(i));
    }
  }
  delete grVIBigFiltered;
  delete grVQBigFiltered;

  // Now need to add noise (if noise terms exist)
  if (!noiseTerms.empty()) {
    std::cout << "Adding noise...\n";
    AddNoise();
  }
}

rad::Signal::Signal(TString filePath, ICavity* cav, LocalOscillator lo,
                    double sRate, std::vector<GaussianNoise> noiseTerms,
                    double tAcq)
    : localOsc(lo), sampleRate(sRate), noiseVec(noiseTerms), cavity(cav) {
  CreateVoltageGraphs();

  // Check if input file opens properly
  SetUpTree(filePath);

  // Set file info
  GetFileInfo();

  // TO DO: figure out what the relevant modes are

  // Calculate mode normalisation
  // This should be for all the relevant modes but for now just do it for the
  // TE111 mode
  const std::complex<double> normTE111p{
      cavity->GetModeNormalisation(ICavity::kTE, 1, 1, 1, true)};
  const std::complex<double> normTE111m{
      cavity->GetModeNormalisation(ICavity::kTE, 1, 1, 1, false)};
  std::cout << "Normalisation +ve, -ve = " << normTE111p << ", " << normTE111p
            << std::endl;

  // Figure out where we're going to generate the signal up to
  // By default, just do the whole electron trajectory file
  if (tAcq < 0) tAcq = fileEndTime;

  double sampleTime{0};    // Second sample time in seconds
  double sample10Time{0};  // First sample time in seconds
  unsigned int sample10Num{0};
  double sample10StepSize{1 / (10 * sRate)};

  // Create just one deque for our
  advancedTimeVec.push_back(std::deque<double>());

  // Create some intermediate level graphs before final sampling
  auto grVIInter = new TGraph();
  auto grVQInter = new TGraph();
  auto grVIBigFiltered = new TGraph();
  auto grVQBigFiltered = new TGraph();

  // Loop through tree entries
  // Initially we are just doing the the higher frequency sampling
  double printTime{0};  // seconds
  const double printInterval{5e-6};
  for (unsigned int iE{0}; iE < inputTree->GetEntries(); iE++) {
    inputTree->GetEntry(iE);
    const double entryTime{time};

    // Check we are still within the acquisition time, stop otherwise
    if (entryTime > tAcq) break;

    if (entryTime >= printTime) {
      std::cout << printTime * 1e6 << " us signal processed...\n";
      printTime += printInterval;
    }

    AddNewCavityTimes(entryTime, TVector3(xPos, yPos, zPos));

    // Do we need to sample now?
    if (entryTime >= sample10Time) {
      // Yes we do
      // First get the retarded time to calculate the fields at
      double tr{GetRetardedTime(sample10Time, 0)};

      // Calculate the mode field amplitudes
      double ei_p{CalcCavityEField(tr, normTE111p).X()};
      double eq_p{ei_p};
      double ei_m{CalcCavityEField(tr, normTE111m).X()};
      double eq_m{ei_m};

      DownmixVoltages(ei_p, eq_p, sample10Time);

      grVIInter->SetPoint(grVIInter->GetN(), sample10Time, ei_p);
      grVQInter->SetPoint(grVQInter->GetN(), sample10Time, eq_p);

      sample10Num++;
      sample10Time = double(sample10Num) * sample10StepSize;

      // We want to filter this signal if it's long enough
      // Pick a good number to do FFTs with (2^n preferably)
      if (grVIInter->GetN() == 32768) {
        auto grVISmallFiltered = BandPassFilter(grVIInter, 0, sRate / 2);
        auto grVQSmallFiltered = BandPassFilter(grVQInter, 0, sRate / 2);
        // Now add these points to the existing graph
        for (int iF{0}; iF < grVISmallFiltered->GetN(); iF++) {
          grVIBigFiltered->SetPoint(grVIBigFiltered->GetN(),
                                    grVISmallFiltered->GetPointX(iF),
                                    grVISmallFiltered->GetPointY(iF));
          grVQBigFiltered->SetPoint(grVQBigFiltered->GetN(),
                                    grVQSmallFiltered->GetPointX(iF),
                                    grVQSmallFiltered->GetPointY(iF));
        }
        delete grVISmallFiltered;
        delete grVQSmallFiltered;
        // Clear the current graphs and start again
        delete grVIInter;
        delete grVQInter;
        grVIInter = new TGraph();
        grVQInter = new TGraph();
      }
    } else {
      continue;
    }
  }

  // Can now safely close the input file
  CloseInputFile();
}

rad::Signal::~Signal() {
  delete grVITime;
  delete grVQTime;
}

void rad::Signal::GetFileInfo() {
  inputTree->GetEntry(0);
  fileStartTime = time;
  inputTree->GetEntry(inputTree->GetEntries() - 1);
  fileEndTime = time;
  filePntsPerTime =
      double(inputTree->GetEntries()) / (fileEndTime - fileStartTime);
}

double rad::Signal::CalcVoltage(double tr, IAntenna* ant) {
  if (tr == -1) {
    // The voltage is from before the signal has reached the antenna
    return 0;
  } else {
    // We actually have to calculate the voltage
    // Start off with a first guess
    int firstGuessTInd{int(round(filePntsPerTime * (tr - fileStartTime)))};

    // Find the appropriate time
    inputTree->GetEntry(firstGuessTInd);
    double firstGuessTime{time};
    unsigned int correctIndex{0};
    if (firstGuessTime == tr) {
      // Easy, no need for interpolation
      TVector3 pos(xPos, yPos, zPos);
      TVector3 vel(xVel, yVel, zVel);
      TVector3 acc(xAcc, yAcc, zAcc);
      ROOT::Math::XYZVector eField{
          CalcEField(ant->GetAntennaPosition(), pos, vel, acc)};
      TVector3 eField2(eField.X(), eField.Y(), eField.Z());
      double voltage{(eField2.Dot(antenna[0]->GetETheta(pos)) +
                      eField2.Dot(antenna[0]->GetEPhi(pos))) *
                     antenna[0]->GetHEff()};
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
      ROOT::Math::XYZVector eField{
          CalcEField(antenna[0]->GetAntennaPosition(), pos, vel, acc)};
      TVector3 eField2(eField.X(), eField.Y(), eField.Z());
      double voltage{(eField2.Dot(antenna[0]->GetETheta(pos)) +
                      eField2.Dot(antenna[0]->GetEPhi(pos))) *
                     antenna[0]->GetHEff()};
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
      ROOT::Math::XYZVector eField{
          CalcEField(antenna[0]->GetAntennaPosition(), pos, vel, acc)};
      TVector3 eField2(eField.X(), eField.Y(), eField.Z());
      double voltage{(eField2.Dot(antenna[0]->GetETheta(pos)) +
                      eField2.Dot(antenna[0]->GetEPhi(pos))) *
                     antenna[0]->GetHEff()};
      voltage /= 2.0;
      vVals.at(iEl) = voltage;
    }

    // Now actually do the cubic interpolation
    double vInterp{CubicInterpolation(timeVals, vVals, tr)};
    return vInterp;
  }
}

TVector3 rad::Signal::CalcCavityEField(double tr, std::complex<double> norm) {
  if (tr == -1) {
    // The voltage is from before the signal has reached the antenna
    return TVector3(0, 0, 0);
  } else {
    // We actually have to calculate the voltage
    // Start off with a first guess
    int firstGuessTInd{int(round(filePntsPerTime * (tr - fileStartTime)))};

    // Find the appropriate time
    inputTree->GetEntry(firstGuessTInd);
    double firstGuessTime{time};
    unsigned int correctIndex{0};
    if (firstGuessTime == tr) {
      // Easy, no need for interpolation
      TVector3 pos(xPos, yPos, zPos);
      ComplexVector3 modeFieldPlus{cavity->GetModalEField(
          pos, ICavity::kTE, norm.real(), 1, 1, 1, true)};
      ComplexVector3 modeFieldMinus{cavity->GetModalEField(
          pos, ICavity::kTE, norm.real(), 1, 1, 1, false)};
      // Calculate the current density
      ComplexVector3 J(xVel, yVel, zVel);
      J *= -TMath::Qe();

      // Now calculate the field amplitudes
      // Assume we're at the resonance
      const double factor{-190};
      std::complex<double> fieldAmpPlus{factor * J.Dot(modeFieldPlus)};
      std::complex<double> fieldAmpMinus{factor * J.Dot(modeFieldMinus)};

      // Now calculate the actual electric field at the probe
      ComplexVector3 probeFieldPlus{
          cavity->GetModalEField(cavity->GetProbePosition(), ICavity::kTE,
                                 norm.real(), 1, 1, 1, true) *
          fieldAmpPlus};
      ComplexVector3 probeFieldMinus{
          cavity->GetModalEField(cavity->GetProbePosition(), ICavity::kTE,
                                 norm.real(), 1, 1, 1, false) *
          fieldAmpMinus};
      return (probeFieldPlus + probeFieldMinus).Real();
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
    } else {
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

    // Now can do some interpolation if we haven't already found the value
    std::vector<double> timeVals(4);
    std::vector<double> ExVals(4);
    std::vector<double> EyVals(4);
    std::vector<double> EzVals(4);
    if (correctIndex == 0) {
      timeVals.at(0) = 0;
      ExVals.at(0) = 0;
      EyVals.at(0) = 0;
      EzVals.at(0) = 0;
    } else {
      inputTree->GetEntry(correctIndex - 1);
      timeVals.at(0) = time;
      TVector3 pos(xPos, yPos, zPos);
      ComplexVector3 modeFieldPlus{cavity->GetModalEField(
          pos, ICavity::kTE, norm.real(), 1, 1, 1, true)};
      ComplexVector3 modeFieldMinus{cavity->GetModalEField(
          pos, ICavity::kTE, norm.real(), 1, 1, 1, false)};
      // Calculate the current density
      ComplexVector3 J(xVel, yVel, zVel);
      J *= -TMath::Qe();

      // Now calculate the field amplitudes
      // Assume we're at the resonance
      const double factor{-190};
      std::complex<double> fieldAmpPlus{factor * J.Dot(modeFieldPlus)};
      std::complex<double> fieldAmpMinus{factor * J.Dot(modeFieldMinus)};
      // Now calculate the actual electric field at the probe
      ComplexVector3 probeFieldPlus{
          cavity->GetModalEField(cavity->GetProbePosition(), ICavity::kTE,
                                 norm.real(), 1, 1, 1, true)};
      probeFieldPlus *= fieldAmpPlus;
      ComplexVector3 probeFieldMinus{
          cavity->GetModalEField(cavity->GetProbePosition(), ICavity::kTE,
                                 norm.real(), 1, 1, 1, false)};
      probeFieldMinus *= fieldAmpMinus;
      TVector3 totalProbeField{(probeFieldPlus + probeFieldMinus).Real()};
      ExVals.at(0) = totalProbeField.X();
      EyVals.at(0) = totalProbeField.Y();
      EzVals.at(0) = totalProbeField.Z();
    }

    // Now add the rest of the elements of the vector
    for (unsigned int iEl{1}; iEl <= 3; iEl++) {
      inputTree->GetEntry(correctIndex + iEl - 1);
      timeVals.at(iEl) = time;
      TVector3 pos(xPos, yPos, zPos);
      ComplexVector3 modeFieldPlus{cavity->GetModalEField(
          pos, ICavity::kTE, norm.real(), 1, 1, 1, true)};
      ComplexVector3 modeFieldMinus{cavity->GetModalEField(
          pos, ICavity::kTE, norm.real(), 1, 1, 1, false)};
      ComplexVector3 J(xVel, yVel, zVel);
      J *= -TMath::Qe();

      // Now calculate the field amplitudes
      // Assume we're at the resonance
      const double factor{-190};
      std::complex<double> fieldAmpPlus{factor * J.Dot(modeFieldPlus)};
      std::complex<double> fieldAmpMinus{factor * J.Dot(modeFieldMinus)};
      // Now calculate the actual electric field at the probe
      ComplexVector3 probeFieldPlus{
          cavity->GetModalEField(cavity->GetProbePosition(), ICavity::kTE,
                                 norm.real(), 1, 1, 1, true)};
      probeFieldPlus *= fieldAmpPlus;
      ComplexVector3 probeFieldMinus{
          cavity->GetModalEField(cavity->GetProbePosition(), ICavity::kTE,
                                 norm.real(), 1, 1, 1, false)};
      probeFieldMinus *= fieldAmpMinus;
      TVector3 totalProbeField{(probeFieldPlus + probeFieldMinus).Real()};
      ExVals.at(iEl) = totalProbeField.X();
      EyVals.at(iEl) = totalProbeField.Y();
      EzVals.at(iEl) = totalProbeField.Z();
    }

    // Now do the cubic interpolation
    double ExInterp{CubicInterpolation(timeVals, ExVals, tr)};
    double EyInterp{CubicInterpolation(timeVals, EyVals, tr)};
    double EzInterp{CubicInterpolation(timeVals, EzVals, tr)};
    return TVector3(ExInterp, EyInterp, EzInterp);
  }
}

void rad::Signal::AddNewTimes(double time, TVector3 ePos) {
  timeVec.push_back(time);

  // Now calculate advanced time for each antenna point
  for (size_t i{0}; i < antenna.size(); i++) {
    double ta{time +
              (ePos - antenna.at(i)->GetAntennaPosition()).Mag() / TMath::C()};
    advancedTimeVec.at(i).push_back(ta);
  }

  // If the deques are getting too long, get rid of the first element
  if (timeVec.size() > 10000) {
    timeVec.pop_front();
    for (size_t i{0}; i < antenna.size(); i++) {
      advancedTimeVec.at(i).pop_front();
    }
  }
}

void rad::Signal::AddNewCavityTimes(double time, TVector3 ePos) {
  timeVec.push_back(time);

  // Calculate the advanced time for cavity probe position
  double ta{time + (ePos - cavity->GetProbePosition()).Mag() / TMath::C()};
  advancedTimeVec.at(0).push_back(ta);

  // If the deques are getting too long, get rid of the first element
  if (timeVec.size() > 10000) {
    timeVec.pop_front();
    advancedTimeVec.at(0).pop_front();
  }
}

double rad::Signal::GetRetardedTime(double ts, unsigned int antInd) {
  unsigned int chosenInd{0};

  // Check if the sample time is before the signal has reached the antenna
  if (ts < advancedTimeVec.at(antInd).at(0)) {
    return -1;
  } else {
    // Find the values which the sample time lives in between
    // Firstly get a good first guess of where to start looking
    unsigned int firstGuess{GetFirstGuessPoint(ts, antInd)};
    // Check which direction to search in
    if (advancedTimeVec.at(antInd).at(firstGuess) < ts) {
      // We are searching upwards
      for (size_t i{firstGuess}; i < advancedTimeVec.at(antInd).size() - 2;
           i++) {
        // Aim to find the advanced time values the ssample point lies between
        if (ts > advancedTimeVec.at(antInd).at(i) &&
            ts < advancedTimeVec.at(antInd).at(i + 1)) {
          chosenInd = i;
          break;
        }
      }
    } else {
      // We are searching downwards
      for (size_t i{firstGuess}; i >= 0; i--) {
        if (ts > advancedTimeVec.at(antInd).at(i) &&
            ts < advancedTimeVec.at(antInd).at(i + 1)) {
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
      advancedTimeVals.at(0) = advancedTimeVec.at(antInd).at(chosenInd - 1);
    }
    timeVals.at(1) = timeVec.at(chosenInd);
    timeVals.at(2) = timeVec.at(chosenInd + 1);
    timeVals.at(3) = timeVec.at(chosenInd + 2);
    advancedTimeVals.at(1) = advancedTimeVec.at(antInd).at(chosenInd);
    advancedTimeVals.at(2) = advancedTimeVec.at(antInd).at(chosenInd + 1);
    advancedTimeVals.at(3) = advancedTimeVec.at(antInd).at(chosenInd + 2);
    return CubicInterpolation(advancedTimeVals, timeVals, ts);
  }
}

void rad::Signal::SetUpTree(TString filePath) {
  // Open the input file
  OpenInputFile(filePath);

  inputTree = (TTree*)inputFile->Get("tree");
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

void rad::Signal::OpenInputFile(TString filePath) {
  // Check if input file opens properly
  inputFile = std::make_unique<TFile>(filePath, "READ");
  if (!inputFile->IsOpen()) {
    std::cout << "Couldn't open file! Exiting.\n";
    exit(1);
  }
}

unsigned int rad::Signal::GetFirstGuessPoint(double ts, unsigned int antInd) {
  const size_t taVecSize{advancedTimeVec.at(antInd).size()};
  const double pntsPerTime{double(taVecSize) /
                           (advancedTimeVec.at(antInd).at(taVecSize - 1) -
                            advancedTimeVec.at(antInd).at(0))};

  size_t firstGuessPnt{
      size_t(round(pntsPerTime * (ts - advancedTimeVec.at(antInd).at(0))))};
  if (firstGuessPnt < 0) firstGuessPnt = 0;
  return firstGuessPnt;
}

void rad::Signal::CloseInputFile() {
  delete inputTree;
  inputFile->Close();
}

void rad::Signal::DownmixVoltages(double& vi, double& vq, double t) {
  vi *= localOsc.GetInPhaseComponent(t);
  vq *= localOsc.GetQuadratureComponent(t);
}

void rad::Signal::AddNoise() {
  // Set up the noise terms
  for (auto& n : noiseVec) {
    n.SetSampleFreq(sampleRate);
    n.SetSigma();
  }

  // Now actually add the noise
  for (int i{0}; i < grVITime->GetN(); i++) {
    double vi{grVITime->GetPointY(i)};
    double vq{grVQTime->GetPointY(i)};
    for (auto& n : noiseVec) {
      vi += n.GetNoiseVoltage(true);
      vq += n.GetNoiseVoltage(true);
    }
    grVITime->SetPointY(i, vi);
    grVQTime->SetPointY(i, vq);
  }
}

TGraph* rad::Signal::GetVIPowerPeriodogram(double loadResistance) {
  TGraph* grOut = MakePowerSpectrumPeriodogram(grVITime);
  setGraphAttr(grOut);
  grOut->SetTitle("V_{I}; Frequency [Hz]; Power [W]");
  ScaleGraph(grOut, 1 / loadResistance);
  return grOut;
}

TGraph* rad::Signal::GetVQPowerPeriodogram(double loadResistance) {
  TGraph* grOut = MakePowerSpectrumPeriodogram(grVQTime);
  setGraphAttr(grOut);
  grOut->SetTitle("V_{Q}; Frequency [Hz]; Power [W]");
  ScaleGraph(grOut, 1 / loadResistance);
  return grOut;
}

void rad::Signal::CreateVoltageGraphs() {
  grVITime = new TGraph();
  grVQTime = new TGraph();
  setGraphAttr(grVITime);
  setGraphAttr(grVQTime);
  grVITime->GetYaxis()->SetTitle("V_{I}");
  grVQTime->GetYaxis()->SetTitle("V_{Q}");
  grVITime->GetXaxis()->SetTitle("Time [s]");
  grVQTime->GetXaxis()->SetTitle("Time [s]");
}