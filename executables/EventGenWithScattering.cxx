/*
  EventGenWithScattering.cxx

  Generate realistic events with some scattering for Heer
*/

#include <getopt.h>
#include <unistd.h>

#include <boost/lexical_cast.hpp>
#include <boost/uuid/uuid.hpp>
#include <boost/uuid/uuid_generators.hpp>
#include <boost/uuid/uuid_io.hpp>
#include <cmath>
#include <ctime>
#include <filesystem>
#include <iostream>
#include <random>
#include <string>

#include "BasicFunctions/BasicFunctions.h"
#include "BasicFunctions/Constants.h"
#include "ElectronDynamics/BorisSolver.h"
#include "ElectronDynamics/QTNMFields.h"
#include "H5Cpp.h"
#include "Scattering/ElasticScatter.h"
#include "Scattering/InelasticScatter.h"
#include "SignalProcessing/LocalOscillator.h"
#include "SignalProcessing/Signal.h"
#include "TFile.h"
#include "TGraph.h"
#include "TMultiGraph.h"
#include "TString.h"
#include "TSystem.h"
#include "TTree.h"
#include "TVector3.h"
#include "Waveguides/CircularWaveguide.h"
#include "Waveguides/Probe.h"
#include "Waveguides/WaveguideMode.h"

using namespace rad;
using std::cout;
using std::endl;

// Class for storing electron information
class ElectronInfo {
 public:
  double startTime;   // seconds
  double startKE;     // eV
  double startF;      // Hertz
  double axialF;      // Hertz
  double pitchAngle;  // degrees
  double zMax;        // metres
  double deltaE;      // Energy loss from scatter, eV
  TVector3 startPos;  // metres
  TVector3 startVel;  // m/s

  // Default constructor
  ElectronInfo()
      : startTime(0),
        startKE(-1),
        startF(-1),
        axialF(-1),
        pitchAngle(DBL_MAX),
        zMax(-DBL_MAX),
        deltaE(-1),
        startPos(0, 0, 0),
        startVel(0, 0, 0) {}

  void Reset() {
    startTime = 0;
    startKE = -1;
    startF = -1;
    axialF = -1;
    pitchAngle = DBL_MAX;
    zMax = -DBL_MAX;
    deltaE = -1;
    startPos = TVector3(0, 0, 0);
    startVel = TVector3(0, 0, 0);
  }
};

double GenerateElectron(TString file, TVector3 pos, TVector3 vel,
                        BaseField *field, double gasDensity,
                        long double stepSize, double maxSimTime,
                        std::vector<ElectronInfo> &ei, BaseField *bField,
                        double heFraction) {
  bool isTrapped{true};
  TFile fT(file, "recreate");
  TTree tree("tree", "tree");
  double time{};
  double xPos{}, yPos{}, zPos{};
  double xVel{}, yVel{}, zVel{};
  double xAcc{}, yAcc{}, zAcc{};
  tree.Branch("time", &time);
  tree.Branch("xPos", &xPos);
  tree.Branch("yPos", &yPos);
  tree.Branch("zPos", &zPos);
  tree.Branch("xVel", &xVel);
  tree.Branch("yVel", &yVel);
  tree.Branch("zVel", &zVel);
  tree.Branch("xAcc", &xAcc);
  tree.Branch("yAcc", &yAcc);
  tree.Branch("zAcc", &zAcc);
  time = 0;
  xPos = pos.X();
  yPos = pos.Y();
  zPos = pos.Z();
  xVel = pos.X();
  yVel = pos.Y();
  zVel = pos.Z();
  const double tau{2 * R_E / (3 * TMath::C())};
  BorisSolver solver(field, -TMath::Qe(), ME, tau);
  TVector3 acc{solver.acc(pos, vel)};
  xAcc = acc.X();
  yAcc = acc.Y();
  zAcc = acc.Z();
  tree.Fill();

  TVector3 p{pos};
  TVector3 v{vel};

  // Set up random number stuff
  std::random_device rd;
  std::mt19937 gen(rd());

  // Set some information about the electron creation
  double gamma{1 / sqrt(1 - pow(v.Mag() / TMath::C(), 2))};
  double ke{(gamma - 1) * ME_EV};
  ElectronInfo eiTemp;
  eiTemp.startTime = time;
  eiTemp.startKE = ke;
  eiTemp.startPos = pos;
  eiTemp.startVel = vel;
  eiTemp.deltaE = -1;

  // Set up some things for axial frequency calculation
  double zCrossTimes[3];
  uint zCrossings{0};
  double oldZ{pos.Z()};

  // Set up some things for start frequency calculation
  const double fieldMeasurementTime{1e-6};  // seconds
  double currentInfoWriteTime{0};
  double bMean{0};
  uint nFieldPoints{0};

  bool haveSavedInfo{false};

  unsigned long nSteps{1};

  while (time <= maxSimTime && isTrapped) {
    if (std::fmod(time, 5e-6) < stepSize) {
      std::cout << "Simulated " << time * 1e6 << " us\n";
    }

    // Calculate the total cross-section
    // We want to see if the scatter is happening in the next time step
    double gamma{1 / sqrt(1 - pow(v.Mag() / TMath::C(), 2))};
    double ke{(gamma - 1) * ME_EV};
    ElasticScatter scatterElasticHe(ke, 2, 4);
    ElasticScatter scatterElasticT(ke, 1, 3);
    InelasticScatter scatterInelasticHe(ke, He);
    InelasticScatter scatterInelasticT(ke, H);
    double elasticXSecHe{scatterElasticHe.GetTotalXSec()};
    double inelasticXSecHe{scatterInelasticHe.GetTotalXSec()};
    double totalXSecHe{(elasticXSecHe + inelasticXSecHe) * heFraction};
    double elasticXSecT{scatterElasticT.GetTotalXSec()};
    double inelasticXSecT{scatterInelasticT.GetTotalXSec()};
    double totalXSecT{(elasticXSecT + inelasticXSecT) * (1 - heFraction)};

    // Using the gas fraction, calculate the total cross-section
    double totalXSec{totalXSecHe + totalXSecT};
    double meanFreeTime{1 / (gasDensity * totalXSec * v.Mag())};
    // Draw from an exponential distribution to see if we scatter
    std::exponential_distribution<double> scatterDistro(1 / meanFreeTime);
    double scatterTime{scatterDistro(gen)};
    if (scatterTime < stepSize) {
      time += scatterTime;
      cout << "Scattering after " << time * 1e6 << " us\n";
      currentInfoWriteTime += scatterTime;

      // Move the particle up to the scattering time
      // Advance the electron
      auto outputStep = solver.advance_step(scatterTime, p, v);
      p = std::get<0>(outputStep);
      v = std::get<1>(outputStep);

      // Reset some variables used for calculating truth info
      currentInfoWriteTime = 0;
      zCrossings = 0;
      oldZ = zPos;
      bMean = 0;
      nFieldPoints = 0;
      eiTemp.Reset();
      haveSavedInfo = false;

      eiTemp.startTime = time;
      eiTemp.startKE = ke;
      eiTemp.startPos = p;

      // We now need to determine what we are scattering on and if the scatter
      // is elastic or inelastic
      double scatterAngle{0};
      std::uniform_real_distribution<double> uni(0, 1);
      if (uni(gen) < totalXSecHe / totalXSec) {
        // Scattering from Helium
        // Now check if this is an elastic or inelastic scatter
        if (uni(gen) < elasticXSecHe / (elasticXSecHe + inelasticXSecHe)) {
          // Scattering from Helium
          cout << "Elastic scattering from Helium\n";
          scatterAngle = scatterElasticHe.GetRandomScatteringAngle();
          double newKE{scatterElasticHe.GetEnergyAfterScatter(scatterAngle)};
          eiTemp.deltaE = ke - newKE;
          ke = newKE;
          v = scatterElasticHe.GetScatteredVector(v, ke, scatterAngle);
        } else {
          cout << "Inelastic scattering from Helium\n";
          double W{scatterInelasticHe.GetRandomW()};
          eiTemp.deltaE = ke - W;
          scatterAngle = scatterInelasticHe.GetRandomTheta(W);
          ke = W;
          v = scatterInelasticHe.GetScatteredVector(v, ke, scatterAngle);
        }
      } else {
        // Scattering from Tritium
        // Check if this an elastic or inelastic scatter
        if (uni(gen) < elasticXSecT / (elasticXSecT + inelasticXSecT)) {
          cout << "Elastic scattering from Tritium\n";
          scatterAngle = scatterElasticT.GetRandomScatteringAngle();
          double newKE{scatterElasticT.GetEnergyAfterScatter(scatterAngle)};
          eiTemp.deltaE = ke - newKE;
          ke = newKE;
          v = scatterElasticT.GetScatteredVector(v, ke, scatterAngle);
        } else {
          cout << "Inelastic scattering from Tritium\n";
          double W{scatterInelasticT.GetRandomW()};
          eiTemp.deltaE = ke - W;
          scatterAngle = scatterInelasticT.GetRandomTheta(W);
          ke = W;
          v = scatterInelasticT.GetScatteredVector(v, ke, scatterAngle);
        }
      }

      acc = solver.acc(p, v);
      xPos = p.X();
      yPos = p.Y();
      zPos = p.Z();
      xVel = v.X();
      yVel = v.Y();
      zVel = v.Z();
      xAcc = acc.X();
      yAcc = acc.Y();
      zAcc = acc.Z();
      tree.Fill();

      double pitchAngleDeg{abs(atan(v.Perp() / v.Z())) * 180 / M_PI};
      cout << "Scattering angle = " << scatterAngle * 180 / M_PI
           << " degrees\tNew KE (delta E) = " << ke / 1e3 << " keV ("
           << eiTemp.deltaE << " eV)\tNew pitch angle = " << pitchAngleDeg
           << " degrees\n";
      eiTemp.startVel = v;

    } else {
      time = (long double)nSteps * stepSize;
      currentInfoWriteTime += stepSize;

      // Advance the electron
      auto outputStep = solver.advance_step(stepSize, p, v);
      p = std::get<0>(outputStep);
      v = std::get<1>(outputStep);

      // Check if the electron has escaped
      if (abs(p.Z()) > 7.5e-2) {
        std::cout << "Electron escaped after " << time * 1e6 << " us\n";
        isTrapped = false;

        // If we haven't saved the information yet for the previous scattering
        // then do so now
        if (!haveSavedInfo) {
          // Calculate start frequency
          gamma = 1 / sqrt(1 - pow(v.Mag() / TMath::C(), 2));
          ke = (gamma - 1) * ME_EV;
          bMean /= double(nFieldPoints);
          double startF{CalcCyclotronFreq(ke, bMean)};
          eiTemp.startF = startF;
          // We may or may not be able to do our axial frequency calculation
          cout << "Electron escaped after only " << zCrossings
               << " z crossings.\n";
          if (zCrossings == 2) {
            double axFreq{1 / (2 * (zCrossTimes[1] - zCrossTimes[0]))};
            eiTemp.axialF = axFreq;
          }

          ei.push_back(eiTemp);
          haveSavedInfo = true;
        }

        break;
      }

      // Write to tree
      acc = solver.acc(p, v);
      xPos = p.X();
      yPos = p.Y();
      zPos = p.Z();
      xVel = v.X();
      yVel = v.Y();
      zVel = v.Z();
      xAcc = acc.X();
      yAcc = acc.Y();
      zAcc = acc.Z();
      tree.Fill();

      if (zCrossings < 3) {
        if ((oldZ > 0 && zPos < 0) || (oldZ < 0 && zPos > 0)) {
          zCrossTimes[zCrossings] = time;
          zCrossings++;
        }
        oldZ = zPos;
      }

      // We don't want to calculate all the time
      if (currentInfoWriteTime < fieldMeasurementTime) {
        bMean += bField->evaluate_field_magnitude(p);
        nFieldPoints++;

        // Calculate the instantaneous pitch angle
        double pitchAngleDeg{abs(atan(v.Perp() / v.Z())) * 180 / M_PI};
        if (pitchAngleDeg < eiTemp.pitchAngle)
          eiTemp.pitchAngle = pitchAngleDeg;

        // Update the maximum z position if necessary
        if (zPos > eiTemp.zMax) eiTemp.zMax = zPos;
      } else if (!haveSavedInfo && time >= fieldMeasurementTime) {
        // Now calculate the axial frequency
        bMean /= double(nFieldPoints);
        double startF{CalcCyclotronFreq(ke, bMean)};
        double axFreq{1 / (zCrossTimes[2] - zCrossTimes[0])};
        cout << "f_c = " << startF / 1e9 << " GHz\tf_a = " << axFreq / 1e6
             << " MHz\tz_max = " << eiTemp.zMax * 100
             << " cm\tPitch angle = " << eiTemp.pitchAngle << " degrees\n";
        eiTemp.axialF = axFreq;
        eiTemp.startF = startF;
        ei.push_back(eiTemp);
        haveSavedInfo = true;
      }
      oldZ = zPos;

      nSteps++;
    }
  }

  fT.cd();
  tree.Write("", TObject::kOverwrite);
  fT.Close();
  return time;
}

/// @brief Function for generating a random file string
/// @return A unique string which can be used for file names
std::string make_uuid() {
  return boost::lexical_cast<std::string>((boost::uuids::random_generator())());
}

int main(int argc, char *argv[]) {
  // Start a clock
  const clock_t startTotalClock = clock();

  int opt{};
  std::string outputStemStr{" "};  // Directory in which to store files
  unsigned int nEvents{50};        // Number of events to generate
  double gasDensity{1e18};         // m^-3
  double maxSimTime{10e-3};        // seconds
  bool useBathtub{false};          // Use a bathtub field
  double wgRadius{6.0e-3};         // metres
  double heFraction{0};            // Fraction of helium in the gas
  double minPitchAngle{
      0};                  // Minimum initial pitch angle we simulate in degrees
  double trapDepth{4e-3};  // Tesla, depth of the magnetic trap

  // Long options
  static struct option longOptions[] = {
      {"output", required_argument, nullptr, 'o'},
      {"nevents", required_argument, nullptr, 'n'},
      {"density", required_argument, nullptr, 'd'},
      {"max_time", required_argument, nullptr, 't'},
      {"wg_radius", required_argument, nullptr, 'r'},
      {"he_fraction", required_argument, nullptr, 'f'},
      {"min_pitch_angle", required_argument, nullptr, 'p'},
      {"trap_depth", required_argument, nullptr, 0},
      {"bathtub", no_argument, nullptr, 'b'},
      {"help", no_argument, nullptr, 'h'},
      {nullptr, 0, nullptr, 0}};

  int optionIndex{0};
  while ((opt = getopt_long(argc, argv, ":o:n:d:t:r:f:p:bh", longOptions,
                            &optionIndex)) != -1) {
    switch (opt) {
      case 0:  // Handle the long-only option
        if (std::string(longOptions[optionIndex].name) == "trap_depth") {
          trapDepth = boost::lexical_cast<double>(optarg);
        }
        break;
      case 'o':
        outputStemStr = optarg;
        break;
      case 'n':
        nEvents = boost::lexical_cast<unsigned int>(optarg);
        break;
      case 'd':
        gasDensity = boost::lexical_cast<double>(optarg);
        break;
      case 't':
        maxSimTime = boost::lexical_cast<double>(optarg);
        break;
      case 'r':
        wgRadius = boost::lexical_cast<double>(optarg);
        break;
      case 'f':
        heFraction = boost::lexical_cast<double>(optarg);
        break;
      case 'p':
        minPitchAngle = boost::lexical_cast<double>(optarg);
        break;
      case 'b':
        useBathtub = true;
        break;
      case 'h':
        cout << "Usage: " << argv[0]
             << "[--output|-o] [--nevents|-n] [--density|-d] "
                "[--max_time|-t] [--wg_radius|-r] [--he_fraction|-f]"
                "[--min_pitch_angle|-p] [--bathtub|-b] [--trap_depth]"
                "[--help|-h]"
             << endl;
        return 0;
      case ':':
        std::cerr << "Option -" << static_cast<char>(optopt)
                  << " requires an argument.";
        return 1;
      default:
        std::cerr << "Unknown option " << static_cast<char>(optopt);
        return 1;
    }
  }

  cout << "Writing output files to " << outputStemStr << endl;
  cout << "Attempting to generate " << nEvents << " events" << endl;
  cout << "Maximum simulation time = " << maxSimTime * 1e6 << " us" << endl;
  cout << "Gas density = " << gasDensity << " m^-3" << endl;
  cout << "Trap depth = " << trapDepth << " T" << endl;
  cout << "Helium fraction = " << heFraction << endl;
  cout << "Minimum initial pitch angle = " << minPitchAngle << " degrees"
       << endl;

  // Check if the chosen output directory exists
  bool dirExists{std::filesystem::is_directory(outputStemStr)};
  // If it doesn't exist then attempt to create the directory
  if (!dirExists) {
    bool dirCreated{std::filesystem::create_directories(outputStemStr)};
    if (!dirCreated) {
      cout << "Failed to create the output directory. Exiting.\n";
      exit(1);
    }
  }

  TString outputDir{outputStemStr};

  // Magnitude of background field
  double bkgField{1.0};  // Tesla
  // Define the field depending on the trap type selected
  const double rCoil{20e-3};                        // metres
  const double iCoil{2 * trapDepth * rCoil / MU0};  // Amps
  const double trapLength{0.1};                     // metres, bathtub only

  BaseField *field{nullptr};
  if (useBathtub) {
    // Define the bathtub field
    field = new BathtubField(rCoil, iCoil, -trapLength / 2, trapLength / 2,
                             TVector3(0, 0, bkgField));
  } else {
    // Define the harmonic field
    bkgField += trapDepth;
    field = new HarmonicField(rCoil, iCoil, bkgField);
  }
  TVector3 centralField{field->evaluate_field_at_point(TVector3(0, 0, 0))};
  const double centralFieldMag{centralField.Mag()};
  const double endpointKE{18.575e3};  // eV
  const double centralCycFreq{CalcCyclotronFreq(endpointKE, centralFieldMag)};

  // Do about 15 time steps per cyclotron orbit
  const double simStepSize{1 / (10 * centralCycFreq)};  // seconds

  // Now define the waveguide that we want to collect our signals with
  const double wgLength{20e-2};             // metres
  TVector3 probePos1{0, 0, wgLength / 2};   // Place probe at end of guide
  TVector3 probePos2{0, 0, -wgLength / 2};  // Place probe at end of guide
  Probe probe1_zplus(probePos1, WaveguideMode(1, 1, kTE), true);
  Probe probe2_zplus(probePos1, WaveguideMode(1, 1, kTE), false);
  Probe probe1_zminus(probePos2, WaveguideMode(1, 1, kTE), true);
  Probe probe2_zminus(probePos2, WaveguideMode(1, 1, kTE), false);
  auto wg = new CircularWaveguide(wgRadius, wgLength);

  Probe probeArr[2] = {probe1_zplus, probe2_zplus};
  const char *probeNames[2] = {"1", "2"};
  const unsigned int nProbes{2};

  // Signal processing stuff
  const double sampleRate{1e9};  // Hz
  double extraOffset = (useBathtub) ? 50e6 : 0;
  const double loFreq{centralCycFreq - sampleRate / 4 + extraOffset};  // Hz
  // Define local oscillator
  LocalOscillator lo(2 * M_PI * loFreq);

  // Create random number stuff
  std::random_device rd;
  std::mt19937 gen(rd());

  // Define some simulation parameters
  uint nGenerated{0};
  const double rhoGenMax{0.5 * wgRadius};  // metres
  const double energyWindow{1000};         // eV
  const double zLimit{wgLength / 8};
  while (nGenerated < nEvents) {
    const clock_t startEventClock{clock()};

    // Create a uniform distribution
    std::uniform_real_distribution<double> uni1(0, 1);
    // Generate the initial electron position
    double zGen{0};
    const double rhoGen{rhoGenMax * sqrt(uni1(gen))};
    const double phiPosGen{uni1(gen) * 2 * M_PI};
    TVector3 initPos(rhoGen * cos(phiPosGen), rhoGen * sin(phiPosGen), zGen);
    double xStart{initPos.X()};
    double yStart{initPos.Y()};
    double zStart{initPos.Z()};

    // Generate the initial direction
    const double phiVelGen{uni1(gen) * 2 * M_PI};
    const double thetaVelGen{acos(2 * uni1(gen) - 1)};

    // Generate the initial energy
    const double initialKE{(endpointKE - energyWindow) +
                           energyWindow * uni1(gen)};
    const double initialSpeed{GetSpeedFromKE(initialKE, ME)};
    TVector3 initVel(initialSpeed * sin(thetaVelGen) * cos(phiVelGen),
                     initialSpeed * sin(thetaVelGen) * sin(phiVelGen),
                     initialSpeed * cos(thetaVelGen));
    // Calculate pitch angle
    const double pitchAngleRad{abs(atan(initVel.Perp() / initVel.Z()))};
    const double pitchAngleDeg{pitchAngleRad * 180 / M_PI};

    if (pitchAngleDeg < minPitchAngle) continue;

    // Now generate the trajectory
    std::string trackFileExt{make_uuid()};
    TString trackFile{outputDir + Form("/track_%s.root", trackFileExt.data())};

    std::vector<ElectronInfo> eiVec{};
    double simTime{GenerateElectron(trackFile, initPos, initVel, field,
                                    gasDensity, simStepSize, maxSimTime, eiVec,
                                    field, heFraction)};

    if (simTime > 1e-6) {
      // Print out the information for the different scattering events
      for (auto &iInfo : eiVec) {
        cout << "Start time = " << iInfo.startTime * 1e6
             << " us\tStart KE = " << iInfo.startKE
             << " eV\tStart frequency = " << iInfo.startF / 1e9
             << " GHz\tAxial frequency = " << iInfo.axialF / 1e6
             << " MHz\tMax z = " << iInfo.zMax * 100
             << " cm\tPitch angle = " << iInfo.pitchAngle << " degrees\n";
        cout << "Start position = (" << iInfo.startPos.X() << ", "
             << iInfo.startPos.Y() << ", " << iInfo.startPos.Z()
             << ")\tStart velocity = (" << iInfo.startVel.X() << ", "
             << iInfo.startVel.Y() << ", " << iInfo.startVel.Z() << ")\n\n";
      }

      const clock_t endPropClock{clock()};
      const double propTime{(endPropClock - startEventClock) /
                            (double)CLOCKS_PER_SEC};
      cout << "Event " << nGenerated << " took " << propTime
           << " seconds to propagate electron.\n";

      // Set up one output file per event. This prevents us losing all the
      // data if one event fails.
      std::string outfileExt{make_uuid()};
      std::string FILE_NAME{outputStemStr + "/out_" + outfileExt + ".h5"};
      // Open the file
      auto file = new H5::H5File(FILE_NAME, H5F_ACC_TRUNC);
      // Create a group in the file
      std::string GROUP_NAME{"/Data"};
      auto group = new H5::Group(file->createGroup(GROUP_NAME));
      // Create property list for the dataset and set up fill values
      double fillValue{0};
      H5::DSetCreatPropList plist;
      plist.setFillValue(H5::PredType::NATIVE_DOUBLE, &fillValue);

      // Now generate the signals for the probes one at a time
      for (unsigned int iProbe{0}; iProbe < nProbes; iProbe++) {
        // Start a clock for the signal generation
        const clock_t startSignalClock{clock()};

        Signal signal(trackFile, wg, lo, sampleRate, probeArr[iProbe]);
        // Get the output voltage
        auto grV{signal.GetVITimeDomain()};

        // Create dataspace for dataset
        const unsigned int DSPACE_DIM = grV->GetN();
        cout << "Time series is " << DSPACE_DIM << " entries long\n";
        const unsigned int DSPACE_RANK{1};
        hsize_t dim[] = {DSPACE_DIM};
        auto dspace = new H5::DataSpace(DSPACE_RANK, dim);

        // Create dataset and write it into the file
        const std::string DATASET1_NAME(GROUP_NAME + "/signal" +
                                        probeNames[iProbe]);
        auto dataset1 = new H5::DataSet(file->createDataSet(
            DATASET1_NAME, H5::PredType::NATIVE_DOUBLE, *dspace, plist));

        // Define attributes
        // First create the attribute for the time step
        H5::DataSpace timeStepSpace(H5S_SCALAR);
        H5::Attribute timeStepAttr1{dataset1->createAttribute(
            "Time step [seconds]", H5::PredType::NATIVE_DOUBLE, timeStepSpace)};
        double timeStep{1.0 / sampleRate};
        timeStepAttr1.write(H5::PredType::NATIVE_DOUBLE, &timeStep);
        // Local oscillator frequency
        H5::DataSpace loFreqSpc(H5S_SCALAR);
        H5::Attribute loFreqAttr1{dataset1->createAttribute(
            "LO frequency [Hertz]", H5::PredType::NATIVE_DOUBLE, loFreqSpc)};
        loFreqAttr1.write(H5::PredType::NATIVE_DOUBLE, &loFreq);
        // Write the metadata for the trap config
        H5::DataSpace rCoilSpc(H5S_SCALAR);
        H5::Attribute rCoilAttr1{dataset1->createAttribute(
            "r_coil [metres]", H5::PredType::NATIVE_DOUBLE, rCoilSpc)};
        rCoilAttr1.write(H5::PredType::NATIVE_DOUBLE, &rCoil);

        H5::DataSpace iCoilSpc(H5S_SCALAR);
        H5::Attribute iCoilAttr1{dataset1->createAttribute(
            "i_coil [Amps]", H5::PredType::NATIVE_DOUBLE, iCoilSpc)};
        iCoilAttr1.write(H5::PredType::NATIVE_DOUBLE, &iCoil);

        H5::DataSpace bkgFieldSpc(H5S_SCALAR);
        H5::Attribute bkgFieldAttr1{dataset1->createAttribute(
            "B_bkg [Tesla]", H5::PredType::NATIVE_DOUBLE, bkgFieldSpc)};
        bkgFieldAttr1.write(H5::PredType::NATIVE_DOUBLE, &bkgField);

        // Write the metadata for the waveguide dimenstions as well
        H5::DataSpace rWgSpc(H5S_SCALAR);
        H5::Attribute rWgAttr1{dataset1->createAttribute(
            "r_wg [metres]", H5::PredType::NATIVE_DOUBLE, rWgSpc)};
        rWgAttr1.write(H5::PredType::NATIVE_DOUBLE, &wgRadius);

        // Metadata for gas density
        H5::DataSpace gasDensitySpc(H5S_SCALAR);
        H5::Attribute gasDensityAttr1{dataset1->createAttribute(
            "Gas density [m^-3]", H5::PredType::NATIVE_DOUBLE, gasDensitySpc)};
        gasDensityAttr1.write(H5::PredType::NATIVE_DOUBLE, &gasDensity);

        // Metadata for the helium fraction
        H5::DataSpace heFractionSpc(H5S_SCALAR);
        H5::Attribute heFractionAttr1{dataset1->createAttribute(
            "Helium fraction", H5::PredType::NATIVE_DOUBLE, heFractionSpc)};
        heFractionAttr1.write(H5::PredType::NATIVE_DOUBLE, &heFraction);

        // Now create the attributes for each scattering event
        const unsigned int SCATTERSPACE_DIM = eiVec.size();
        hsize_t dimScatter[] = {SCATTERSPACE_DIM};
        H5::DataSpace scatterSpace(1, dimScatter);
        H5::Attribute timeAttr1{dataset1->createAttribute(
            "Time [seconds]", H5::PredType::NATIVE_DOUBLE, scatterSpace)};
        double timeBuffer[SCATTERSPACE_DIM];

        H5::Attribute startEAttr1{dataset1->createAttribute(
            "Start E [eV]", H5::PredType::NATIVE_DOUBLE, scatterSpace)};
        double startEBuffer[SCATTERSPACE_DIM];

        H5::Attribute startFAttr1{dataset1->createAttribute(
            "Start F [Hertz]", H5::PredType::NATIVE_DOUBLE, scatterSpace)};
        double startFBuffer[SCATTERSPACE_DIM];

        H5::Attribute startFDMAttr1{dataset1->createAttribute(
            "Start F, downmixed [Hertz]", H5::PredType::NATIVE_DOUBLE,
            scatterSpace)};
        double startFDMBuffer[SCATTERSPACE_DIM];

        H5::Attribute axialFAttr1{dataset1->createAttribute(
            "Axial F [Hertz]", H5::PredType::NATIVE_DOUBLE, scatterSpace)};
        double axialFBuffer[SCATTERSPACE_DIM];

        H5::Attribute zMaxAttr1{dataset1->createAttribute(
            "z_max [metres]", H5::PredType::NATIVE_DOUBLE, scatterSpace)};
        double zMaxBuffer[SCATTERSPACE_DIM];

        H5::Attribute pitchAngleAttr1{dataset1->createAttribute(
            "Pitch angle [degrees]", H5::PredType::NATIVE_DOUBLE,
            scatterSpace)};
        double pitchAngleBuffer[SCATTERSPACE_DIM];

        H5::Attribute xpStartAttr1{dataset1->createAttribute(
            "x_start [metres]", H5::PredType::NATIVE_DOUBLE, scatterSpace)};
        double xpStartBuffer[SCATTERSPACE_DIM];

        H5::Attribute ypStartAttr1{dataset1->createAttribute(
            "y_start [metres]", H5::PredType::NATIVE_DOUBLE, scatterSpace)};
        double ypStartBuffer[SCATTERSPACE_DIM];

        H5::Attribute zpStartAttr1{dataset1->createAttribute(
            "z_start [metres]", H5::PredType::NATIVE_DOUBLE, scatterSpace)};
        double zpStartBuffer[SCATTERSPACE_DIM];

        H5::Attribute deltaEAttr1{dataset1->createAttribute(
            "Delta E [eV]", H5::PredType::NATIVE_DOUBLE, scatterSpace)};
        double deltaEBuffer[SCATTERSPACE_DIM];

        // Calculate the impedance of the waveguide and write it as an
        // attribute
        const double zWg{wg->GetModeImpedance(WaveguideMode(1, 1, kTE),
                                              2 * M_PI * eiVec[0].startF)};
        H5::DataSpace zWgSpc(H5S_SCALAR);
        H5::Attribute zWgAttr1{dataset1->createAttribute(
            "Waveguide impedance [Ohms]", H5::PredType::NATIVE_DOUBLE, zWgSpc)};
        zWgAttr1.write(H5::PredType::NATIVE_DOUBLE, &zWg);

        // Write the load impedance of 50 Ohms
        const double loadImp{50};  // Ohms
        H5::DataSpace loadImpSpc(H5S_SCALAR);
        H5::Attribute loadImpAttr{dataset1->createAttribute(
            "Load impedance [Ohms]", H5::PredType::NATIVE_DOUBLE, loadImpSpc)};
        loadImpAttr.write(H5::PredType::NATIVE_DOUBLE, &loadImp);

        // Stop the signal generation clock
        const clock_t endSignalClock{clock()};
        const double signalGenTime{(endSignalClock - startSignalClock) /
                                   (double)CLOCKS_PER_SEC};
        // Write the generation times to file
        H5::DataSpace ePropTimeSpc(H5S_SCALAR);
        H5::Attribute ePropAttr{dataset1->createAttribute(
            "ePropTime", H5::PredType::NATIVE_DOUBLE, ePropTimeSpc)};
        ePropAttr.write(H5::PredType::NATIVE_DOUBLE, &propTime);

        H5::DataSpace signalGenTimeSpc(H5S_SCALAR);
        H5::Attribute signalGenAttr{dataset1->createAttribute(
            "signalGenTime", H5::PredType::NATIVE_DOUBLE, signalGenTimeSpc)};
        signalGenAttr.write(H5::PredType::NATIVE_DOUBLE, &signalGenTime);

        for (uint i{0}; i < eiVec.size(); i++) {
          timeBuffer[i] = eiVec[i].startTime;
          startFBuffer[i] = eiVec[i].startF;
          startFDMBuffer[i] = eiVec[i].startF - loFreq;
          axialFBuffer[i] = eiVec[i].axialF;
          startEBuffer[i] = eiVec[i].startKE;
          zMaxBuffer[i] = eiVec[i].zMax;
          pitchAngleBuffer[i] = eiVec[i].pitchAngle;
          xpStartBuffer[i] = eiVec[i].startPos.X();
          ypStartBuffer[i] = eiVec[i].startPos.Y();
          zpStartBuffer[i] = eiVec[i].startPos.Z();
          deltaEBuffer[i] = eiVec[i].deltaE;
        }
        timeAttr1.write(H5::PredType::NATIVE_DOUBLE, &timeBuffer);
        startEAttr1.write(H5::PredType::NATIVE_DOUBLE, &startEBuffer);
        startFAttr1.write(H5::PredType::NATIVE_DOUBLE, &startFBuffer);
        startFDMAttr1.write(H5::PredType::NATIVE_DOUBLE, &startFDMBuffer);
        axialFAttr1.write(H5::PredType::NATIVE_DOUBLE, &axialFBuffer);
        zMaxAttr1.write(H5::PredType::NATIVE_DOUBLE, &zMaxBuffer);
        pitchAngleAttr1.write(H5::PredType::NATIVE_DOUBLE, &pitchAngleBuffer);
        xpStartAttr1.write(H5::PredType::NATIVE_DOUBLE, &xpStartBuffer);
        ypStartAttr1.write(H5::PredType::NATIVE_DOUBLE, &ypStartBuffer);
        zpStartAttr1.write(H5::PredType::NATIVE_DOUBLE, &zpStartBuffer);
        deltaEAttr1.write(H5::PredType::NATIVE_DOUBLE, &deltaEBuffer);

        // Create the buffer for writing in the voltage data
        double bufferIn[DSPACE_DIM];
        for (uint i{0}; i < grV->GetN(); i++) {
          bufferIn[i] = grV->GetPointY(i);
        }
        dataset1->write(bufferIn, H5::PredType::NATIVE_DOUBLE, *dspace);

        // Close the dataset and dataspace
        delete dataset1;
        delete dspace;
      }

      // Close the group
      delete group;
      // Close the HDF5 file
      delete file;

      const clock_t endEventClock{clock()};
      const double eventTime{(endEventClock - startEventClock) /
                             (double)CLOCKS_PER_SEC};
      cout << "Event " << nGenerated << " took " << eventTime
           << " seconds to generate.\n";
      nGenerated++;
    }
    // Delete the track file
    cout << "Deleting track file\n\n";
    gSystem->Exec("rm -f " + trackFile);
  }
  cout << "Finished generating all events.\n";

  delete field;
  delete wg;
  return 0;
}
