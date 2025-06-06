/*
  MLSampleGeneration.cxx

  This source code make an executable which generates an ensemble of events for
  a given pitch angle range and radial distribution and energy. Outputs are
  written to a HDF5 file.
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
#include "ElectronDynamics/TrajectoryGen.h"
#include "H5Cpp.h"
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

/// @brief Function for generating a random file string
/// @return A unique string which can be used for file names
std::string make_uuid() {
  return boost::lexical_cast<std::string>((boost::uuids::random_generator())());
}

double GetTrueCycFreq(const TString &trackFile, BaseField *field) {
  // Open the track file
  TFile file(trackFile);
  TTree *tree = (TTree *)file.Get("tree");
  double time{0};
  double xPos{0}, yPos{0}, zPos{0};
  double xVel{0}, yVel{0}, zVel{0};
  double bMean{0.0};
  tree->SetBranchAddress("time", &time);
  tree->SetBranchAddress("xPos", &xPos);
  tree->SetBranchAddress("yPos", &yPos);
  tree->SetBranchAddress("zPos", &zPos);
  tree->SetBranchAddress("xVel", &xVel);
  tree->SetBranchAddress("yVel", &yVel);
  tree->SetBranchAddress("zVel", &zVel);
  tree->GetEntry(0);

  unsigned int nPoints{1};
  bMean += field->evaluate_field_magnitude(TVector3(xPos, yPos, zPos));

  const double maxTime{1e-6};
  while (time < maxTime) {
    tree->GetEntry(nPoints);
    bMean += field->evaluate_field_magnitude(TVector3(xPos, yPos, zPos));
    nPoints++;
  }
  bMean /= double(nPoints);

  tree->GetEntry(0);
  TVector3 initVel(xVel, yVel, zVel);
  const double gamma{1 / sqrt(1 - pow(initVel.Mag() / TMath::C(), 2))};
  const double EInit{(gamma - 1) * ME_EV};

  delete tree;
  file.Close();

  return CalcCyclotronFreq(EInit, bMean);
}

double GetTrueAxialFreq(const TString &trackFile) {
  TFile file(trackFile);
  TTree *tree = (TTree *)file.Get("tree");
  double time{0};
  double zPos{0};
  tree->SetBranchAddress("time", &time);
  tree->SetBranchAddress("zPos", &zPos);
  tree->GetEntry(0);

  double zCrossTimes[3];
  unsigned int zCrossings{0};
  double oldZ{zPos};

  int iEntry{1};
  while (zCrossings < 3) {
    tree->GetEntry(iEntry);
    if ((oldZ > 0 && zPos < 0) || (oldZ < 0 && zPos > 0)) {
      zCrossTimes[zCrossings] = time;
      zCrossings++;
    }
    oldZ = zPos;
    iEntry++;
  }

  delete tree;
  file.Close();
  return 1.0 / (zCrossTimes[2] - zCrossTimes[0]);
}

int main(int argc, char *argv[]) {
  // Start a clock
  const clock_t startTotalClock = clock();

  int opt{};
  std::string outputStemStr{" "};       // Directory in which to store files
  unsigned int nEvents{50};             // Number of events to generate
  double maxSimTime{10e-6};             // seconds
  double minPitchAngle{89.5};           // degrees
  double energyWindow{0};               // eV
  bool bathtubTrap{false};              // Use a bathtub trap
  double bkgField{1.0};                 // Tesla
  double wgRadius{5.0e-3};              // metres
  double xMax{2e-3};                    // Max start radius in metres
  const double endpointEnergy{18.6e3};  // eV

  while ((opt = getopt(argc, argv, ":o:n:e:p:t:f:r:x:bh")) != -1) {
    switch (opt) {
      case 'o':
        outputStemStr = optarg;
        break;
      case 'n':
        nEvents = boost::lexical_cast<unsigned int>(optarg);
        break;
      case 'e':
        energyWindow = boost::lexical_cast<double>(optarg);
        break;
      case 'p':
        minPitchAngle = boost::lexical_cast<double>(optarg);
        break;
      case 't':
        maxSimTime = boost::lexical_cast<double>(optarg);
        break;
      case 'f':
        bkgField = boost::lexical_cast<double>(optarg);
        break;
      case 'r':
        wgRadius = boost::lexical_cast<double>(optarg);
        break;
      case 'x':
        xMax = boost::lexical_cast<double>(optarg);
        break;
      case 'b':
        bathtubTrap = true;
        break;
      case ':':
        cout << "Option -" << static_cast<char>(optopt)
             << " requires an argument." << endl;
        return 1;
      case '?':
        cout << "Unknown option: " << static_cast<char>(optopt) << endl;
        return 1;
      case 'h':
        cout << "Usage: " << argv[0]
             << " [-o output directory] [-n number of events] [-e energy] [-p "
                "min pitch angle] [-t simulation time] [-f background field] "
                "[-r waveguide radius] [-x max electron start radius] [-b "
                "bathtub trap]"
             << endl;
        return 1;
    }
  }

  cout << "Writing output files to " << outputStemStr << "\n";
  cout << "Attempting to generate " << nEvents << " events.\n";
  if (energyWindow > 0) {
    cout << "Generating electron energies between "
         << endpointEnergy - energyWindow << " and " << endpointEnergy << " eV"
         << endl;
  } else {
    cout << "Kinetic energy of electrons = " << endpointEnergy << " eV" << endl;
  }

  cout << "Minimum pitch angle of electrons = " << minPitchAngle
       << " degrees\n";
  cout << "Background field = " << bkgField << " Tesla\n";
  cout << "Simulation time = " << maxSimTime * 1e6 << " us\n";
  cout << "Waveguide radius = " << wgRadius * 1e3 << " mm\n";
  cout << "Maximumum start radius = " << xMax * 1e3 << " mm\n";
  if (bathtubTrap) {
    cout << "Using bathtub trap\n";
  } else {
    cout << "Using harmonic trap\n";
  }

  // Check if waveguide radius is less than the start radius
  if (wgRadius < xMax) {
    cout
        << "Waveguide radius must be greater than the start radius. Exiting.\n";
    exit(1);
  }

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

  // Create a directory for the temporary track files
  std::string trackDir{outputStemStr + "/tracks"};
  bool trackDirExists{std::filesystem::is_directory(trackDir)};
  if (!trackDirExists) {
    bool trackDirCreated{std::filesystem::create_directories(trackDir)};
    if (!trackDirCreated) {
      cout << "Failed to create the track directory. Exiting.\n";
      exit(1);
    }
  }

  TString outputDir{outputStemStr};

  // Define the field depending on the trap type selected
  const double rCoil{2e-2};                         // m
  const double deltaTheta{3.9866802 * M_PI / 180};  // radians
  const double trapDepth{bkgField *
                         (1 / pow(cos(deltaTheta), 2) - 1)};  // Tesla
  cout << "Trap depth = " << trapDepth << " Tesla" << endl;
  const double iCoil{2 * trapDepth * rCoil / MU0};  // Amps
  const double trapLength{0.1};                     // metres, bathtub only
  BaseField *field{nullptr};
  if (bathtubTrap) {
    // Create the field for the bathtub trap
    field = new BathtubField(rCoil, iCoil, -trapLength / 2, trapLength / 2,
                             TVector3(0, 0, bkgField));
  } else {
    // Create the field for the harmonic trap
    bkgField += trapDepth;
    field = new HarmonicField(rCoil, iCoil, bkgField);
  }

  // Calculate the rough cyclotron frequency
  const double centralCycFreq{CalcCyclotronFreq(endpointEnergy, bkgField)};

  // Set number of steps per cyclotron orbit
  const double simStepSize{1 / (10 * centralCycFreq)};

  // Now define the waveguide that we want to collect our signal with
  const double wgLength{20e-2};            // metres
  TVector3 probePos1{0, 0, wgLength / 2};  // Place probe at end of guide

  auto wg = new CircularWaveguide(wgRadius, wgLength);

  // Signal processing things
  const double sampleRate{1e9};  // Hz
  double extraOffset = (bathtubTrap) ? 50e6 : 0;
  const double loFreq{centralCycFreq - sampleRate / 4 + extraOffset};  // Hz
  // Define local oscillator
  LocalOscillator lo(2 * M_PI * loFreq);

  // Create random number stuff
  std::random_device rd;
  std::mt19937 gen(rd());

  // Draw pitch angle from a uniform distribution
  std::uniform_real_distribution<double> pitchDist(minPitchAngle, 90.0);

  // Start generating events
  for (unsigned int iEvent{0}; iEvent < nEvents; iEvent++) {
    const clock_t startEventClock{clock()};
    // Create a uniform distribution
    std::uniform_real_distribution<double> uni1(0, 1);
    // Generation parameters
    const double energy{endpointEnergy - energyWindow +
                        energyWindow * uni1(gen)};
    const double initialSpeed{GetSpeedFromKE(energy, ME)};

    cout << "Generating event " << iEvent + 1 << " of " << nEvents << ": "
         << "Energy = " << energy / 1e3 << " keV\n";

    const double pitchAngle{pitchDist(gen)};
    const double pitchAngleRadians{pitchAngle * M_PI / 180};

    // Generate a random position on a disk
    const double rhoGen{xMax * sqrt(uni1(gen))};
    const double phiPosGen{uni1(gen) * 2 * M_PI};
    TVector3 initPos(rhoGen * cos(phiPosGen), rhoGen * sin(phiPosGen), 0);
    double xStart{initPos.X()};
    double yStart{initPos.Y()};
    double zStart{initPos.Z()};

    // Generate the initial velocity. The pitch angle is now fixed but vary the
    // cyclotron phase
    const double cycPhase{uni1(gen) * 2 * M_PI};
    // Randomly choose positive or negative z velocity
    double zComp{1};
    if (uni1(gen) < 0.5) {
      zComp = -1;
    }
    TVector3 initVel(sin(pitchAngleRadians) * cos(cycPhase),
                     sin(pitchAngleRadians) * sin(cycPhase),
                     zComp * cos(pitchAngleRadians));
    initVel *= initialSpeed;

    // Now generate the trajectory
    std::string trackFileExt{make_uuid()};
    TString trackFile{trackDir + Form("/track_%s.root", trackFileExt.data())};
    ElectronTrajectoryGen traj(trackFile, field, initPos, initVel, simStepSize,
                               maxSimTime);

    double trueCycFreq{GetTrueCycFreq(trackFile, field)};
    cout << "True cyclotron frequency = " << trueCycFreq / 1e6 << " MHz\n";
    double trueAxialFreq{GetTrueAxialFreq(trackFile)};
    cout << "True axial frequency = " << trueAxialFreq / 1e6 << " MHz\n";

    // Set up one output file per event. This prevents us losing all the data
    // if one event fails.
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

    // Now generate the signals. One for each polarisation
    Probe probePlus(probePos1, WaveguideMode(1, 1, kTE), true);
    Probe probeMinus(probePos1, WaveguideMode(1, 1, kTE), false);
    const size_t nProbes{2};
    Probe probeArr[nProbes] = {probePlus, probeMinus};
    for (size_t iProbe{0}; iProbe < nProbes; iProbe++) {
      Signal sig(trackFile, wg, lo, sampleRate, probeArr[iProbe]);
      auto grV{sig.GetVITimeDomain()};

      // Create dataspace for dataset
      const unsigned int DSPACE_DIM = grV->GetN();
      cout << "Time series is " << DSPACE_DIM << " entries long\n";
      const unsigned int DSPACE_RANK{1};
      hsize_t dim[] = {DSPACE_DIM};
      auto dspace = new H5::DataSpace(DSPACE_RANK, dim);

      // Create dataset and write it into the file
      const std::string DATASET_NAME1(GROUP_NAME + "/signal" +
                                      std::to_string(iProbe + 1));
      auto dataset1 = new H5::DataSet(file->createDataSet(
          DATASET_NAME1, H5::PredType::NATIVE_DOUBLE, *dspace, plist));

      // Define attributes for datasets
      // First create the attribute for the time step
      H5::DataSpace timeStepSpace(H5S_SCALAR);
      H5::Attribute timeStepAttr{dataset1->createAttribute(
          "Time step [seconds]", H5::PredType::NATIVE_DOUBLE, timeStepSpace)};
      double timeStep{1.0 / sampleRate};
      timeStepAttr.write(H5::PredType::NATIVE_DOUBLE, &timeStep);
      // Local oscillator frequency
      H5::DataSpace loFreqSpc(H5S_SCALAR);
      H5::Attribute loFreqAttr{dataset1->createAttribute(
          "LO frequency [Hertz]", H5::PredType::NATIVE_DOUBLE, loFreqSpc)};
      loFreqAttr.write(H5::PredType::NATIVE_DOUBLE, &loFreq);

      // Pitch angle information
      H5::DataSpace paSpace(H5S_SCALAR);
      H5::Attribute paAttr{dataset1->createAttribute(
          "Pitch angle [degrees]", H5::PredType::NATIVE_DOUBLE, paSpace)};
      paAttr.write(H5::PredType::NATIVE_DOUBLE, &pitchAngle);

      // Pitch angle information
      H5::DataSpace axFSpace(H5S_SCALAR);
      H5::Attribute axFAttr{dataset1->createAttribute(
          "Axial F [Hertz]", H5::PredType::NATIVE_DOUBLE, axFSpace)};
      axFAttr.write(H5::PredType::NATIVE_DOUBLE, &trueAxialFreq);

      // Energy information
      H5::DataSpace energySpace(H5S_SCALAR);
      H5::Attribute energyAttr{dataset1->createAttribute(
          "Energy [eV]", H5::PredType::NATIVE_DOUBLE, energySpace)};
      energyAttr.write(H5::PredType::NATIVE_DOUBLE, &energy);

      // Starting position information as 3D vector
      H5::DataSpace posSpace(H5S_SIMPLE);
      hsize_t posDim[] = {3};
      posSpace.setExtentSimple(1, posDim);
      H5::Attribute posAttr{dataset1->createAttribute(
          "Starting position [metres]", H5::PredType::NATIVE_DOUBLE, posSpace)};
      double posData[3] = {xStart, yStart, zStart};
      posAttr.write(H5::PredType::NATIVE_DOUBLE, posData);

      // Starting velocity information as 3D vector
      H5::DataSpace velSpace(H5S_SIMPLE);
      hsize_t velDim[] = {3};
      velSpace.setExtentSimple(1, velDim);
      H5::Attribute velAttr{
          dataset1->createAttribute("Starting velocity [metres/second]",
                                    H5::PredType::NATIVE_DOUBLE, velSpace)};
      double velData[3] = {initVel.X(), initVel.Y(), initVel.Z()};
      velAttr.write(H5::PredType::NATIVE_DOUBLE, velData);

      // Write the true cyclotron frequency
      H5::DataSpace trueCycFreqSpc(H5S_SCALAR);
      H5::Attribute trueCycFreqAttr{dataset1->createAttribute(
          "Cyclotron frequency [Hertz]", H5::PredType::NATIVE_DOUBLE,
          trueCycFreqSpc)};
      trueCycFreqAttr.write(H5::PredType::NATIVE_DOUBLE, &trueCycFreq);

      // Write the downmixed cyclotron frequency
      H5::DataSpace downmixCycFreqSpc(H5S_SCALAR);
      H5::Attribute downmixCycFreqAttr{dataset1->createAttribute(
          "Downmixed cyclotron frequency [Hertz]", H5::PredType::NATIVE_DOUBLE,
          downmixCycFreqSpc)};
      double downmixCycFreq{trueCycFreq - loFreq};
      downmixCycFreqAttr.write(H5::PredType::NATIVE_DOUBLE, &downmixCycFreq);

      // Calculate the impedance of the waveguide and write as an attribute
      WaveguideMode te11Mode(1, 1, kTE);
      const double zWg{wg->GetModeImpedance(te11Mode, 2 * M_PI * trueCycFreq)};
      cout << "Mode impedance = " << zWg << " Ohms\n";
      H5::DataSpace zWgSpc(H5S_SCALAR);
      H5::Attribute zWgAttr{dataset1->createAttribute(
          "Waveguide impedance [Ohms]", H5::PredType::NATIVE_DOUBLE, zWgSpc)};
      zWgAttr.write(H5::PredType::NATIVE_DOUBLE, &zWg);

      // Write the load impedance of 50 Ohms
      const double loadImp{50};  // Ohms
      H5::DataSpace loadImpSpc(H5S_SCALAR);
      H5::Attribute loadImpAttr{dataset1->createAttribute(
          "Load impedance [Ohms]", H5::PredType::NATIVE_DOUBLE, loadImpSpc)};
      loadImpAttr.write(H5::PredType::NATIVE_DOUBLE, &loadImp);

      // Write the metadata for the trap config
      H5::DataSpace rCoilSpc(H5S_SCALAR);
      H5::Attribute rCoilAttr{dataset1->createAttribute(
          "r_coil [metres]", H5::PredType::NATIVE_DOUBLE, rCoilSpc)};
      rCoilAttr.write(H5::PredType::NATIVE_DOUBLE, &rCoil);

      H5::DataSpace iCoilSpc(H5S_SCALAR);
      H5::Attribute iCoilAttr{dataset1->createAttribute(
          "i_coil [Amps]", H5::PredType::NATIVE_DOUBLE, iCoilSpc)};
      iCoilAttr.write(H5::PredType::NATIVE_DOUBLE, &iCoil);

      H5::DataSpace bkgFieldSpc(H5S_SCALAR);
      H5::Attribute bkgFieldAttr{dataset1->createAttribute(
          "B_bkg [Tesla]", H5::PredType::NATIVE_DOUBLE, bkgFieldSpc)};
      bkgFieldAttr.write(H5::PredType::NATIVE_DOUBLE, &bkgField);

      // Write the metadata for the waveguide dimenstions as well
      H5::DataSpace rWgSpc(H5S_SCALAR);
      H5::Attribute rWgAttr{dataset1->createAttribute(
          "r_wg [metres]", H5::PredType::NATIVE_DOUBLE, rWgSpc)};
      rWgAttr.write(H5::PredType::NATIVE_DOUBLE, &wgRadius);

      if (bathtubTrap) {
        H5::DataSpace trapLengthSpc(H5S_SCALAR);
        H5::Attribute trapLengthAttr{dataset1->createAttribute(
            "Trap length [metres]", H5::PredType::NATIVE_DOUBLE,
            trapLengthSpc)};
        trapLengthAttr.write(H5::PredType::NATIVE_DOUBLE, &trapLength);
      }

      // Create the buffer for writing in the data
      double bufferIn[DSPACE_DIM];
      for (unsigned int i{0}; i < DSPACE_DIM; i++) {
        bufferIn[i] = grV->GetPointY(i);
      }
      dataset1->write(bufferIn, H5::PredType::NATIVE_DOUBLE, *dspace);

      // Close the dataset and dataspace
      delete dataset1;
      delete dspace;
    }

    // Delete the track file
    gSystem->Exec("rm -f " + trackFile);

    // Close the group
    delete group;
    // Close the HDF5 file
    delete file;

    const clock_t endEventClock = clock();
    cout << "Event generation took "
         << double(endEventClock - startEventClock) / CLOCKS_PER_SEC
         << " s\n\n";
  }
}
