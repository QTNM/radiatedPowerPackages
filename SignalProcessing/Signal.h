/*
  Signal.h

  Rewritten 01/09/23 by S. Jones
*/

#ifndef SIGNAL_H
#define SIGNAL_H

#include <deque>
#include <memory>
#include <vector>

#include "Antennas/IAntenna.h"
#include "BasicFunctions/ButterworthFilter.h"
#include "SignalProcessing/LocalOscillator.h"
#include "SignalProcessing/NoiseFunc.h"
#include "TFile.h"
#include "TGraph.h"
#include "TTree.h"
#include "Waveguides/ICavity.h"
#include "Waveguides/IWaveguide.h"
#include "Waveguides/Probe.h"
#include "Waveguides/WaveguideMode.h"

namespace rad {
class Signal {
 public:
  /// @brief Parametrised constructor
  /// @param trajectoryFilePath String to electron trajectory file
  /// @param ant Pointer to antenna
  /// @param lo Local oscillator
  /// @param sRate Sample rate in Hertz
  /// @param noiseTerms Vector of noise terms
  /// @param tAcq Acquisition time for signal in seconds
  Signal(TString trajectoryFilePath, IAntenna* ant, LocalOscillator lo,
         double sRate, std::vector<GaussianNoise> noiseTerms = {},
         double tAcq = -1);

  /// @brief Parametrised constructor for multiple antennas
  /// @param trajectoryFilePath String to electron trajectory file
  /// @param ant Pointer to antenna
  /// @param lo Local oscillator
  /// @param sRate Sample rate in Hertz
  /// @param noiseTerms Vector of noise terms
  /// @param tAcq Acquisition time for signal in seconds
  Signal(TString trajectoryFilePath, std::vector<IAntenna*> ant,
         LocalOscillator lo, double sRate,
         std::vector<GaussianNoise> noiseTerms = {}, double tAcq = -1);

  /// @brief Parametrised constructor using cavity as RF collection device
  /// @param filePath String to electron trajectory file
  /// @param cav Cavity used for RF collection
  /// @param lo Local oscillator for downmixing
  /// @param sRate Sample rate in Hertz
  /// @param noiseTerms Vector of noise terms
  /// @param tAcq Acquisition time for signal in seconds
  Signal(TString filePath, ICavity* cav, LocalOscillator lo, double sRate,
         std::vector<GaussianNoise> noiseTerms = {}, double tAcq = -1);

  /// @brief Parametrised constructor using waveguide as RF collection device
  /// @param filePath String to electron trajectory file
  /// @param wg Pointer to waveguide
  /// @param lo Local oscillator
  /// @param sRate Sample rate in Hertz
  /// @param probe Probe for signal readout
  /// @param noiseTerms Vector of noise to add
  /// @param tAcq Max acqusition time in seconds
  Signal(TString filePath, IWaveguide* wg, LocalOscillator lo, double sRate,
         Probe probe, std::vector<GaussianNoise> noiseTerms = {},
         double tAcq = -1);

  /// Destructor
  ~Signal();

  /// @brief Getter function for in-phase voltage component. Assumes 50 Ohm load
  /// @return Time domain voltage graph
  TGraph* GetVITimeDomain() { return grVITime; }

  /// @brief Getter function for quadrature voltage component. Assumes 50 Ohm
  /// load
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

  // Simulation step size
  long double simStepSize;

  // Noise terms
  std::vector<GaussianNoise> noiseVec;

  TGraph* grVITime = 0;  // In phase component
  TGraph* grVQTime = 0;  // Quadrature component

  std::deque<long double> timeVec;
  std::vector<std::deque<long double>>
      advancedTimeVec;  // One deque per antenna

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
  std::vector<IAntenna*> antenna;

  // Pointer to cavity if necessary
  ICavity* cavity;

  // Pointer to waveguide if necessary
  IWaveguide* waveguide;

  // Probe for waveguide
  Probe pr;

  /// @brief Function to add new times to vectors
  /// @param time New time from file in seconds
  /// @param ePos Electron position vector in metres
  void AddNewTimes(long double time, TVector3 ePos);

  /// @brief Function to add new times for vector
  /// @param time New time from file in seconds
  /// @param ePos Electron position vector in metres
  /// @param probe Probe to use for signal readout
  /// @param freq Frequency of signal in Hertz
  void AddNewCavWgTimes(long double time, TVector3 ePos, Probe& probe,
                        double freq);

  /// @brief Calculate the retarded time
  /// @param ts Sample time in seconds
  /// @param antInd Index of antenna
  /// @return Relevant retarded time in seconds
  long double GetRetardedTime(long double ts, unsigned int antInd);

  /// @brief Get a guess for the index to start at
  /// @param ts Sample time in seconds
  /// @return Index of guess
  int GetFirstGuessPoint(long double ts, unsigned int antInd);

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
  /// @param ant Pointer to chosen antenna
  /// @return Voltage in volts
  double CalcVoltage(long double tr, IAntenna* ant);

  /// @brief Calculate electric field at cavity probe position
  /// @param tr Retarded time in seconds
  /// @param norm Field normalisation
  /// @return Electric field at point in V/m
  TVector3 CalcCavityEField(double tr, std::complex<double> norm);

  /// @brief Calculate electric field at waveguide probe position
  /// @param tr Retarded time in seconds
  /// @param mode Waveguide mode to calculate
  /// @param norm Field normalisation
  /// @param omega Angular frequency of cyclotron motion
  /// @return Electric field at point in V/m
  TVector3 CalcWaveguideEField(double tr, WaveguideMode mode, double norm,
                               double omega, bool state);

  /// @brief Calculate electric field at waveguide probe position
  /// @param tr Retarded time in seconds
  /// @param mode Waveguide mode to calculate
  /// @param omega Angular frequency of cyclotron motion
  /// @return Waveguide mode amplitude
  double CalcWgAmp(double tr, WaveguideMode mode, double omega);

  /// @brief Function for downmixing voltages
  /// @param vi In phase voltage component
  /// @param vq Quadrature voltage component
  /// @param t Time in seconds
  void DownmixVoltages(double& vi, double& vq, double t);

  /// @brief Sets some key parameters about the input file
  void GetFileInfo();

  /// @brief Adds noise to the signals
  void AddNoise();

  /// @brief Sets up voltage graphs for class
  void CreateVoltageGraphs();

  /// @brief Calculates the initial electron frequency from trajectory info
  /// @return Initial cyclotron frequency in Hertz
  double CalcInitialFreq();
};

}  // namespace rad

#endif