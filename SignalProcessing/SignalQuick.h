/*
  QuickSignal.h
*/

#ifndef SIGNAL_QUICK_H
#define SIGNAL_QUICK_H

#include <memory>
#include <vector>

#include "Antennas/IAntenna.h"
#include "SignalProcessing/LocalOscillator.h"
#include "SignalProcessing/NoiseFunc.h"
#include "TGraph.h"
#include "TTree.h"

namespace rad {
class SignalQuick {
 public:
  SignalQuick(TString trajectoryFilePath, std::shared_ptr<IAntenna> ant,
              LocalOscillator lo, double sRate,
              std::vector<GaussianNoise> noiseTerms = {});

  TGraph GetVITimeDomain() { return grVITime; }

  TGraph GetVQTimeDomain() { return grVQTime; }

 private:
  LocalOscillator localOsc;

  TGraph grVITime;  // In phase component
  TGraph grVQTime;  // Quadrature component

  std::vector<double> timeVec;
  std::vector<double> advancedTimeVec;

  TVector3 antennaPos;

  // Pointer to input file
  std::unique_ptr<TFile> inputFile = 0;
  // Pointer to input tree from file
  TTree *inputTree = 0;
  // Variables for input tree
  double time{};
  double xPos{}, yPos{}, zPos{};
  double xVel{}, yVel{}, zVel{};
  double xAcc{}, yAcc{}, zAcc{};

  // Pointer to the antenna
  std::shared_ptr<IAntenna> antenna = 0;

  /// @brief Function to add new times to vectors
  /// @param time New time from file in seconds
  /// @param ePos Electron position vector in metres
  void AddNewTimes(double time, TVector3 ePos);

  /// @brief Calculate the retarded time
  /// @param ts Sample time in seconds
  /// @return Relevant retarded time in seconds
  double GetRetardedTime(double ts);

  /// @brief Get a guess for the index to start at
  /// @param ts Sample time in seconds
  /// @return Index of guess
  unsigned int GetFirstGuessPoint(double ts);

  /// @brief Sets up the tree to be read in
  /// @param filePath Path to electron trajectory file
  void SetUpTree(TString filePath);

  /// @brief Safely open the input
  /// @param filePath Path to electron trajectory file
  void OpenInputFile(TString filePath);

  /// @brief Safely closes the input file
  void CloseInputFile();

  /// @brief Calculate the voltage at a given time
  /// @param tr Retarded time at which to calculate the voltage [seconds]
  /// @return Voltage in volts
  double CalcVoltage(double tr);
};

}  // namespace rad

#endif