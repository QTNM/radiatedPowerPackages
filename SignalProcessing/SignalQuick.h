/*
  QuickSignal.h
*/

#ifndef SIGNAL_QUICK_H
#define SIGNAL_QUICK_H

#include <deque>
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
  /// @brief Parametrised constructor
  /// @param trajectoryFilePath
  /// @param ant
  /// @param lo
  /// @param sRate
  /// @param noiseTerms Vector of noise terms
  SignalQuick(TString trajectoryFilePath, IAntenna* ant, LocalOscillator lo,
              double sRate, std::vector<GaussianNoise> noiseTerms = {});

  /// Destructor
  ~SignalQuick();

  /// @brief Getter function for in-phase voltage component
  /// @return Time domain voltage graph
  TGraph* GetVITimeDomain() { return grVITime; }

  /// @brief Getter function for quadrature voltage component
  /// @return Time domain voltage graph
  TGraph* GetVQTimeDomain() { return grVQTime; }

  /// @brief Returns the power spectrum of the in-phase voltage component after
  /// all signal processing Power spectrum in this case is the periodogram
  /// @param loadResistance The load resistance used for the power calculation
  /// @return Pointer to power spectrum
  TGraph* GetVIPowerPeriodogram(double loadResistance);

  /// @brief Returns the power spectrum of the quadrature voltage component
  /// after all signal processing Power spectrum in this case is the periodogram
  /// @param loadResistance The load resistance used for the power calculation
  /// @return Pointer to power spectrum
  TGraph* GetVQPowerPeriodogram(double loadResistance);

 private:
  // Oscillator for the downmixing
  LocalOscillator localOsc;

  // Sample rate
  double sampleRate;

  // Noise terms
  std::vector<GaussianNoise> noiseVec;

  TGraph* grVITime = 0;  // In phase component
  TGraph* grVQTime = 0;  // Quadrature component

  std::deque<double> timeVec;
  std::deque<double> advancedTimeVec;

  TVector3 antennaPos;

  // Input file details
  double fileStartTime{};
  double fileEndTime{};
  double filePntsPerTime{};

  // Pointer to input file
  std::unique_ptr<TFile> inputFile = 0;
  // Pointer to input tree from file
  TTree* inputTree = 0;
  // Variables for input tree
  double time{};
  double xPos{}, yPos{}, zPos{};
  double xVel{}, yVel{}, zVel{};
  double xAcc{}, yAcc{}, zAcc{};

  // Pointer to the antenna
  IAntenna* antenna = 0;

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

  /// @brief Function for downmixing voltages
  /// @param vi In phase voltage component
  /// @param vq Quadrature voltage component
  /// @param t Time in seconds
  void DownmixVoltages(double& vi, double& vq, double t);

  /// @brief Sets some key parameters about the input file
  void GetFileInfo();

  /// @brief Adds noise to the signals
  void AddNoise();
};

}  // namespace rad

#endif