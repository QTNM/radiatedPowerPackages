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
  advancedTimeVec.push_back(std::deque<long double>());

  // Create a deque for filtered and unfiltered info
  // Need this for time domain filter
  std::deque<double> unfilteredVi;
  std::deque<double> filteredVi;
  // Copies for quadrature components
  std::deque<double> unfilteredVq;
  std::deque<double> filteredVq;

  // Define the Butterworth filter
  // Want to get rid of frequencies above half the sample rate
  ButterworthFilter filter(6, sRate / 2, 10 * sRate);

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
      long double tr{GetRetardedTime(sample10Time, 0)};

      double vi{CalcVoltage(tr, antenna[0])};
      double vq{vi};

      DownmixVoltages(vi, vq, sample10Time);

      // Now we want to filter this value and add it to the deques
      double viFiltered{filter.FilterValues(unfilteredVi, filteredVi, vi)};
      double vqFiltered{filter.FilterValues(unfilteredVq, filteredVq, vq)};
      // Only need to add 1 in 10 of these to the final signal graph
      if (sample10Num % 10 == 0) {
        grVITime->SetPoint(grVITime->GetN(), sample10Time, viFiltered);
        grVQTime->SetPoint(grVQTime->GetN(), sample10Time, vqFiltered);
      }

      sample10Num++;
      sample10Time = double(sample10Num) * sample10StepSize;
    } else {
      continue;
    }
  }

  // Can now safely close the input file
  CloseInputFile();

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

  double sampleTime{0};                       // Second sample time in seconds
  long double sample10Time{0};                // First sample time in seconds
  unsigned int sample10Num{0};                // First sample number
  double sample10StepSize{1 / (10 * sRate)};  // Initial sampling step size

  // Create number of deques equal to number of antennas
  for (size_t i{0}; i < antenna.size(); i++) {
    advancedTimeVec.push_back(std::deque<long double>());
  }

  // Create a deque for filtered and unfiltered info
  // Need this for time domain filter
  std::deque<double> unfilteredVi;
  std::deque<double> filteredVi;
  // Copies for quadrature components
  std::deque<double> unfilteredVq;
  std::deque<double> filteredVq;

  // Define the Butterworth filter
  // Want to get rid of frequencies above half the sample rate
  ButterworthFilter filter(6, sRate / 2, 10 * sRate);

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
        long double tr{GetRetardedTime(sample10Time, iAnt)};
        double v{CalcVoltage(tr, antenna[iAnt])};
        vi += v;
        vq += v;
      }
      DownmixVoltages(vi, vq, sample10Time);

      // Now we want to filter this value and add it to the deques
      double viFiltered{filter.FilterValues(unfilteredVi, filteredVi, vi)};
      double vqFiltered{filter.FilterValues(unfilteredVq, filteredVq, vq)};
      // Only need to add 1 in 10 of these to the final signal graph
      if (sample10Num % 10 == 0) {
        grVITime->SetPoint(grVITime->GetN(), sample10Time, viFiltered);
        grVQTime->SetPoint(grVQTime->GetN(), sample10Time, vqFiltered);
      }

      sample10Num++;
      sample10Time = double(sample10Num) * sample10StepSize;
    } else {
      continue;
    }
  }

  // Can now safely close the input file
  CloseInputFile();

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

  const double omega{CalcInitialFreq() * TMath::TwoPi()};

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

  // Create just one deque for our probe position
  // To do: allow multiple probes
  advancedTimeVec.push_back(std::deque<long double>());

  // Create a deque for filtered and unfiltered info
  // Need this for time domain filter
  std::deque<double> unfilteredVi;
  std::deque<double> filteredVi;
  // Copies for quadrature components
  std::deque<double> unfilteredVq;
  std::deque<double> filteredVq;

  // Define the Butterworth filter
  // Want to get rid of frequencies above half the sample rate
  ButterworthFilter filter(6, sRate / 2, 10 * sRate);

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

    AddNewCavWgTimes(entryTime, TVector3(xPos, yPos, zPos), pr,
                     omega / (2 * M_PI));

    // Do we need to sample now?
    if (entryTime >= sample10Time) {
      // Yes we do
      // First get the retarded time to calculate the fields at
      long double tr{GetRetardedTime(sample10Time, 0)};

      // Calculate the mode field amplitudes
      double ei_p{CalcCavityEField(tr, normTE111p).X()};
      double eq_p{ei_p};
      double ei_m{CalcCavityEField(tr, normTE111m).X()};
      double eq_m{ei_m};

      DownmixVoltages(ei_p, eq_p, sample10Time);

      // Now we want to filter this value and add it to the deques
      double viFiltered{filter.FilterValues(unfilteredVi, filteredVi, ei_p)};
      double vqFiltered{filter.FilterValues(unfilteredVq, filteredVq, eq_p)};
      // We only need to save 1 out of every 10 of these values to the final
      // signal graphs
      if (sample10Num % 10 == 0) {
        grVITime->SetPoint(grVITime->GetN(), sample10Time, viFiltered);
        grVQTime->SetPoint(grVQTime->GetN(), sample10Time, vqFiltered);
      }

      sample10Num++;
      sample10Time = double(sample10Num) * sample10StepSize;
    } else {
      continue;
    }
  }

  // Can now safely close the input file
  CloseInputFile();

  // Now need to add noise (if noise terms exist)
  if (!noiseTerms.empty()) {
    std::cout << "Adding noise...\n";
    AddNoise();
  }
}

rad::Signal::Signal(TString filePath, IWaveguide* wg, LocalOscillator lo,
                    double sRate, Probe probe,
                    std::vector<GaussianNoise> noiseTerms, double tAcq)
    : localOsc(lo),
      sampleRate(sRate),
      noiseVec(noiseTerms),
      waveguide(wg),
      pr(probe) {
  CreateVoltageGraphs();
  // Check if input file opens properly
  SetUpTree(filePath);
  // Set file info
  GetFileInfo();

  // First things first, we need to figure out a rough frequency our electron is
  // radiating at
  const double omega{CalcInitialFreq() * TMath::TwoPi()};
  std::cout << "Calculated frequency = " << omega * 1e-9 / TMath::TwoPi()
            << " GHz\n";

  // Now need to check if the mode our probe is reading out propagates
  WaveguideMode wm{pr.GetMode()};
  if (!wg->ModePropagates(wm, omega / (2 * M_PI))) {
    std::cerr << "Waveguide mode " << wm
              << "does not propagate. These results cannot be trusted\n ";
    return;
  }

  // Calculate mode integral and resulting normalisation
  const unsigned int nSurfPnts{100};
  double modeNorm{0};
  double integral{0};
  if (waveguide->MultiplePolarisations()) {
    // Check if there are actually multiple polarisations (i.e. they have
    // equal non-zero integrals
    integral = waveguide->GetEFieldIntegral(wm, omega, 1, nSurfPnts,
                                            pr.GetPolarisationState());
    double integralTest2{waveguide->GetEFieldIntegral(
        wm, omega, 1, nSurfPnts, !pr.GetPolarisationState())};
    if ((integral - integralTest2) / integral < 0.01) {
      // We have two polarisations so require an additional factor of two
      integral /= 2;
    } else {
      std::cout << "Only one polarisation found\n";
    }
  } else {
    integral = waveguide->GetEFieldIntegral(wm, omega, 1, nSurfPnts, true);
  }
  if (integral != 0) modeNorm = 1.0 / sqrt(integral);

  // Figure out where we're going to generate the signal up to
  // By default, just do the whole electron trajectory file
  if (tAcq < 0) tAcq = fileEndTime;

  double sampleTime{0};    // Second sample time in seconds
  double sample10Time{0};  // First sample time in seconds
  unsigned int sample10Num{0};
  double sample10StepSize{1 / (10 * sRate)};

  // Create just one deque for our probe position
  // TO DO: Allow for multiple probes
  advancedTimeVec.push_back(std::deque<long double>());

  // Create a deque for filtered and unfiltered info
  // Need this for time domain filter
  std::deque<double> unfilteredVi;
  std::deque<double> filteredVi;
  // Copies for quadrature components
  std::deque<double> unfilteredVq;
  std::deque<double> filteredVq;

  // Define the Butterworth filter
  // Want to get rid of frequencies above half the sample rate
  ButterworthFilter filter(6, sRate / 2, 10 * sRate);

  // Loop through tree entries
  // Initially we are just doing the the higher frequency sampling
  double printTime{0};  // seconds
  // Controls how many seconds are processed before we print
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

    AddNewCavWgTimes(entryTime, TVector3(xPos, yPos, zPos), pr,
                     omega / (2 * M_PI));

    // Do we need to sample now?
    if (entryTime >= sample10Time) {
      // Yes we do
      // First get the retarded time to calculate the fields at
      long double tr{GetRetardedTime(sample10Time, 0)};

      double vi{0};
      double modeImp{waveguide->GetModeImpedance(wm, omega)};
      vi += CalcWgAmp(tr, wm, modeNorm, omega);

      double vq{vi};
      DownmixVoltages(vi, vq, sample10Time);

      // Now we want to filter this value and add it to the deques
      double viFiltered{filter.FilterValues(unfilteredVi, filteredVi, vi)};
      double vqFiltered{filter.FilterValues(unfilteredVq, filteredVq, vq)};
      // We only need to save 1 out of every 10 of these values to the final
      // signal graphs
      if (sample10Num % 10 == 0) {
        grVITime->SetPoint(grVITime->GetN(), sample10Time, viFiltered);
        grVQTime->SetPoint(grVQTime->GetN(), sample10Time, vqFiltered);
      }

      sample10Num++;
      sample10Time = double(sample10Num) * sample10StepSize;
    } else {
      continue;
    }
  }

  // Safely close input file
  CloseInputFile();

  // Now need to add noise (if noise terms exist)
  if (!noiseTerms.empty()) {
    std::cout << "Adding noise...\n";
    AddNoise();
  }
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

  // Now get the simulation step size
  inputTree->GetEntry(1);
  long double time1{time};
  simStepSize = time1 - fileStartTime;
}

double rad::Signal::CalcVoltage(long double tr, IAntenna* ant) {
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
    std::vector<long double> timeVals(4);
    std::vector<long double> vVals(4);
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
    long double vInterp{CubicInterpolation(timeVals, vVals, tr)};
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

TVector3 rad::Signal::CalcWaveguideEField(double tr, WaveguideMode mode,
                                          double norm, double omega,
                                          bool state) {
  if (tr == -1) {
    // The voltage is from before the signal has reached the probe
    // Therefore there is no electric field there at this time
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
      TVector3 pos(xPos, yPos, zPos);
      TVector3 vel(xVel, yVel, zVel);
      // Calculate the field amplitudes
      double ampPlus{
          waveguide->GetFieldAmp(mode, omega, pos, vel, norm, true, true)};
      double ampMinus{
          waveguide->GetFieldAmp(mode, omega, pos, vel, norm, false, true)};
      // Now calculate the actual field at the probe position
      TVector3 probeFieldPlus{
          waveguide->GetModeEField(pos, mode, norm, omega, true)};
      probeFieldPlus *= ampPlus;
      TVector3 probeFieldMinus{
          waveguide->GetModeEField(pos, mode, norm, omega, false)};
      probeFieldMinus *= ampMinus;
      return probeFieldPlus + probeFieldMinus;
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
      TVector3 vel(xVel, yVel, zVel);
      // Calculate the field amplitudes
      double ampPlus{
          waveguide->GetFieldAmp(mode, omega, pos, vel, norm, true, true)};
      double ampMinus{
          waveguide->GetFieldAmp(mode, omega, pos, vel, norm, false, true)};
      // Now calculate the actual field at the probe position
      TVector3 probeFieldPlus{
          waveguide->GetModeEField(pos, mode, norm, omega, true)};
      probeFieldPlus *= ampPlus;
      TVector3 probeFieldMinus{
          waveguide->GetModeEField(pos, mode, norm, omega, false)};
      probeFieldMinus *= ampMinus;
      TVector3 totalProbeField{probeFieldPlus + probeFieldMinus};
      ExVals.at(0) = totalProbeField.X();
      EyVals.at(0) = totalProbeField.Y();
      EzVals.at(0) = totalProbeField.Z();
    }

    // Now add the other elements for the interpolation
    for (unsigned int iEl{1}; iEl <= 3; iEl++) {
      inputTree->GetEntry(correctIndex + iEl - 1);
      timeVals.at(iEl) = time;
      TVector3 pos(xPos, yPos, zPos);
      TVector3 vel(xVel, yVel, zVel);
      // Calculate the field amplitudes
      double ampPlus{
          waveguide->GetFieldAmp(mode, omega, pos, vel, norm, true, true)};
      double ampMinus{
          waveguide->GetFieldAmp(mode, omega, pos, vel, norm, false, true)};

      // Now calculate the actual field at the probe position
      TVector3 probeFieldPlus{
          waveguide->GetModeEField(pos, mode, norm, omega, true)};
      probeFieldPlus *= ampPlus;
      TVector3 probeFieldMinus{
          waveguide->GetModeEField(pos, mode, norm, omega, false)};
      probeFieldMinus *= ampMinus;
      TVector3 totalProbeField{probeFieldPlus + probeFieldMinus};
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

double rad::Signal::CalcWgAmp(double tr, WaveguideMode mode, double norm,
                              double omega) {
  if (tr == -1) {
    // The voltage is from before the signal has reached the probe
    // Therefore there is no electric field there at this time
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
      TVector3 pos(xPos, yPos, zPos);
      TVector3 vel(xVel, yVel, zVel);
      // Calculate the field amplitudes after checking where our probe is in
      // relation to the electron
      if (pr.GetPosition().Z() > pos.Z()) {
        return waveguide->GetFieldAmp(mode, omega, pos, vel, norm,
                                      pr.GetPolarisationState(), true);
      } else {
        return waveguide->GetFieldAmp(mode, omega, pos, vel, norm,
                                      pr.GetPolarisationState(), false);
      }
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
    std::vector<double> ampVals(4);
    if (correctIndex == 0) {
      timeVals.at(0) = 0;
      ampVals.at(0) = 0;
    } else {
      inputTree->GetEntry(correctIndex - 1);
      timeVals.at(0) = time;
      TVector3 pos(xPos, yPos, zPos);
      TVector3 vel(xVel, yVel, zVel);
      // Calculate the field amplitudes after checking where our probe is in
      // relation to the electron
      if (pr.GetPosition().Z() > pos.Z()) {
        ampVals.at(0) = waveguide->GetFieldAmp(mode, omega, pos, vel, norm,
                                               pr.GetPolarisationState(), true);
      } else {
        ampVals.at(0) = waveguide->GetFieldAmp(
            mode, omega, pos, vel, norm, pr.GetPolarisationState(), false);
      }
    }

    // Now add the other elements for the interpolation
    for (unsigned int iEl{1}; iEl <= 3; iEl++) {
      inputTree->GetEntry(correctIndex + iEl - 1);
      timeVals.at(iEl) = time;
      TVector3 pos(xPos, yPos, zPos);
      TVector3 vel(xVel, yVel, zVel);
      // Calculate the field amplitudes
      if (pr.GetPosition().Z() > pos.Z()) {
        ampVals.at(iEl) = waveguide->GetFieldAmp(
            mode, omega, pos, vel, norm, pr.GetPolarisationState(), true);
      } else {
        ampVals.at(iEl) = waveguide->GetFieldAmp(
            mode, omega, pos, vel, norm, pr.GetPolarisationState(), false);
      }
    }

    double ampInterp{CubicInterpolation(timeVals, ampVals, tr)};
    return ampInterp;
  }
}

void rad::Signal::AddNewTimes(long double time, TVector3 ePos) {
  timeVec.push_back(time);

  // Now calculate advanced time for each antenna point
  for (size_t i{0}; i < antenna.size(); i++) {
    long double ta{time + (ePos - antenna.at(i)->GetAntennaPosition()).Mag() /
                              TMath::C()};
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

void rad::Signal::AddNewCavWgTimes(long double time, TVector3 ePos,
                                   Probe& probe, double freq) {
  timeVec.push_back(time);

  // Calculate the advanced time for cavity probe position
  const double v_group{waveguide->GetGroupVelocity(probe.GetMode(), freq)};
  long double ta{time + (ePos - probe.GetPosition()).Mag() / v_group};
  advancedTimeVec.at(0).push_back(ta);

  // If the deques are getting too long, get rid of the first element
  if (timeVec.size() > 10000) {
    timeVec.pop_front();
    advancedTimeVec.at(0).pop_front();
  }
}

long double rad::Signal::GetRetardedTime(long double ts, unsigned int antInd) {
  unsigned int chosenInd{0};

  // Check if the sample time is before the signal has reached the antenna
  if (ts < advancedTimeVec.at(antInd).at(0)) {
    return -1;
  } else {
    // Find the values which the sample time lives in between
    // Firstly get a good first guess of where to start looking
    int firstGuess{GetFirstGuessPoint(ts, antInd)};
    // Some debugging info
    if (firstGuess >= advancedTimeVec.at(antInd).size()) {
      std::cout << "Strange index produced!: firstTime - ts = "
                << (advancedTimeVec.at(antInd).at(0) - ts) << " seconds\n";
    }

    // Check which direction to search in
    if (advancedTimeVec.at(antInd).at(firstGuess) < ts) {
      // We are searching upwards
      for (int i{firstGuess}; i < advancedTimeVec.at(antInd).size() - 2; i++) {
        // Aim to find the advanced time values the ssample point lies between
        if (ts > advancedTimeVec.at(antInd).at(i) &&
            ts < advancedTimeVec.at(antInd).at(i + 1)) {
          chosenInd = i;
          break;
        }
      }
    } else {
      // We are searching downwards
      for (int i{firstGuess}; i >= 0; i--) {
        if (ts > advancedTimeVec.at(antInd).at(i) &&
            ts < advancedTimeVec.at(antInd).at(i + 1)) {
          chosenInd = i;
          break;
        }
      }
    }

    // Now we have the relevant index, we need to interpolate
    // This should be the retarded time
    std::vector<long double> timeVals(4);
    std::vector<long double> advancedTimeVals(4);
    if (chosenInd == 0) {
      // Set these to something sensible
      long double deltaT{timeVec.at(chosenInd + 1) - timeVec.at(chosenInd)};
      long double deltaAT{advancedTimeVals.at(chosenInd + 1) -
                          advancedTimeVals.at(chosenInd)};
      timeVals.at(0) = timeVec.at(chosenInd) - deltaT;
      advancedTimeVals.at(0) =
          advancedTimeVec.at(antInd).at(chosenInd) - deltaAT;
    } else if (chosenInd > advancedTimeVec.at(antInd).size() - 3) {
      chosenInd = advancedTimeVec.at(antInd).size() - 3;
      timeVals.at(0) = timeVec.at(chosenInd - 1);
      advancedTimeVals.at(0) = advancedTimeVec.at(antInd).at(chosenInd - 1);
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

int rad::Signal::GetFirstGuessPoint(long double ts, unsigned int antInd) {
  const size_t taVecSize{advancedTimeVec.at(antInd).size()};
  const long double pntsPerTime{(long double)taVecSize /
                                (advancedTimeVec.at(antInd).at(taVecSize - 1) -
                                 advancedTimeVec.at(antInd).at(0))};

  int firstGuessPnt{
      int(round(pntsPerTime * (ts - advancedTimeVec.at(antInd).at(0))))};
  if (firstGuessPnt < 0)
    return 0;
  else if (firstGuessPnt >= taVecSize)
    return taVecSize - 1;
  else
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

double rad::Signal::CalcInitialFreq() {
  const double maxCalcTime{1e-6};  // seconds
  auto grX = new TGraph();
  for (int i{0}; i < inputTree->GetEntries(); i++) {
    inputTree->GetEntry(i);
    if (time > maxCalcTime) break;

    grX->SetPoint(grX->GetN(), time, xPos);
  }
  auto grXSpec{MakePowerSpectrumPeriodogram(grX)};
  delete grX;

  const double minFreq{1e9};  // Hertz
  double maxPeak{-DBL_MAX};
  double maxPeakFreq{0};
  for (int iX{1}; iX < grXSpec->GetN(); iX++) {
    if (grXSpec->GetPointX(iX) < minFreq) {
      continue;
    } else if (grXSpec->GetPointY(iX) > maxPeak) {
      maxPeak = grXSpec->GetPointY(iX);
      maxPeakFreq = grXSpec->GetPointX(iX);
    }
  }
  delete grXSpec;
  return maxPeakFreq;
}