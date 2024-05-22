/*
  TemplateGeneration.cxx

  Executable allowing for the generation of signal templates
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

using namespace rad;
using std::cout;
using std::endl;

/// @brief Function for generating a random file string
/// @return A unique string which can be used for file names
std::string make_uuid() {
  return boost::lexical_cast<std::string>((boost::uuids::random_generator())());
}

int main(int argc, char* argv[]) {
  // Start a clock
  const clock_t startTotalClock = clock();

  int opt{};
  std::string outputStemStr{" "};  // Directory in which to store files
  double kineticEnergy{18575};     // Electron kinetic energy in eV
  double maxSimTime{10e-6};        // seconds
  while ((opt = getopt(argc, argv, "o:k:t:")) != -1) {
    switch (opt) {
      case 'o':
        outputStemStr = optarg;
        break;
      case 'k':
        kineticEnergy = std::stod(optarg);
        break;
      case 't':
        maxSimTime = boost::lexical_cast<double>(optarg);
        break;
      case ':':
        cout << "Option -" << static_cast<char>(optopt)
             << " requires an argument." << endl;
        return 1;
      case '?':
        cout << "Unknown option: " << static_cast<char>(optopt) << endl;
        return 1;
      default:
        cout << "Usage: " << argv[0]
             << " [-o outputStem] [-k kineticEnergy] [-t templateLength]"
             << endl;
        return 1;
    }
  }

  cout << "Writing output files to " << outputStemStr << "\n";
  cout << "Electron kinetic energy = " << kineticEnergy << " eV\n";
  cout << "Maximum simulation time = " << maxSimTime * 1e6 << " us\n";

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

  // Define the field. Do a harmonic trap.
  const double rCoil{15e-3};                        // metres
  const double trapDepth{3.4e-3};                   // Tesla
  const double iCoil{2 * trapDepth * rCoil / MU0};  // Amps
  // Magnitude of background field
  const double bkgField{0.7};  // Tesla
  auto field = new HarmonicField(rCoil, iCoil, bkgField);
  TVector3 centralField{field->evaluate_field_at_point(TVector3(0, 0, 0))};
  const double centralFieldMag{centralField.Mag()};
  const double endpointKE{18.575e3};  // eV
  const double centralCycFreq{CalcCyclotronFreq(endpointKE, centralFieldMag)};

  // Do about 15 time steps per cyclotron orbit
  const double simStepSize{1 / (10 * centralCycFreq)};  // seconds

  // Now define the waveguide that we want to collect our signal with
  const double wgLength{20e-2};           // metres
  const double wgRadius{7.14e-3};         // metres
  TVector3 probePos{0, 0, wgLength / 2};  // Place probe at end of guide
  auto wg = new CircularWaveguide(wgRadius, wgLength, probePos);

  // Signal processing stuff
  const double sampleRate{1e9};                                  // Hz
  const double loFreq{centralCycFreq - sampleRate / 4 + 150e6};  // Hz
  // Define local oscillator
  LocalOscillator lo(2 * M_PI * loFreq);

  // Create random number stuff
  std::random_device rd;
  std::mt19937 gen(rd());

  // Figure out which points to run
  const double rhoGenMin{0};     // metres
  const double rhoGenMax{6e-3};  // metres
  const uint nRhoPnts{15};
  const double phiGenMin{0};         // radians
  const double phiGenMax{2 * M_PI};  // radians
  const uint nPhiPnts{15};
  const double pitchAngleGenMin{86 * M_PI / 180};  // radians
  const double pitchAngleGenMax{90 * M_PI / 180};  // radians
  const uint nPitchAnglePnts{15};
  const uint nTotalSims{nRhoPnts * nPhiPnts * nPitchAnglePnts -
                        (nPhiPnts - 1) * nPitchAnglePnts};
  uint simCounter{1};

  // Convert kinetic energy to string
  std::string kineticEnergyStr{std::to_string(int(kineticEnergy))};
  // Open an HDF5 file to store the data
  std::string hdf5FileStr{outputStemStr + "/templates_" + kineticEnergyStr +
                          "eV.h5"};
  auto file = new H5::H5File(hdf5FileStr, H5F_ACC_TRUNC);
  // Create a group in the file
  std::string GROUP_NAME{"/templates"};
  auto group = new H5::Group(file->createGroup(GROUP_NAME));
  // Create property list for the dataset and set up fill values
  double fillValue{0};
  H5::DSetCreatPropList plist;
  plist.setFillValue(H5::PredType::NATIVE_DOUBLE, &fillValue);

  // Loop through all of these various points
  for (uint iRho{0}; iRho < nRhoPnts; iRho++) {
    double rho{rhoGenMin +
               double(iRho) * (rhoGenMax - rhoGenMin) / double(nRhoPnts - 1)};
    for (uint iPhi{0}; iPhi < nPhiPnts; iPhi++) {
      double phi{phiGenMin +
                 double(iPhi) * (phiGenMax - phiGenMin) / double(nPhiPnts)};
      if (rho == 0 && phi != 0) continue;

      for (uint iPitch{0}; iPitch < nPitchAnglePnts; iPitch++) {
        const clock_t startEventClock = clock();
        double pitchAngle{pitchAngleGenMin +
                          double(iPitch) *
                              (pitchAngleGenMax - pitchAngleGenMin) /
                              double(nPitchAnglePnts - 1)};

        cout << "Generating electron " << simCounter << " of " << nTotalSims
             << " with (x, y) = (" << rho * cos(phi) << ", " << rho * sin(phi)
             << ")" << " and pitch angle = " << pitchAngle * 180 / M_PI
             << " degrees\n";
        const double zGen{0};  // metres
        TVector3 startPos{rho * cos(phi), rho * sin(phi), zGen};
        double xStart{startPos.X()};
        double yStart{startPos.Y()};
        double initialSpeed{GetSpeedFromKE(kineticEnergy, ME)};
        TVector3 initialVel(sin(pitchAngle), 0, cos(pitchAngle));
        initialVel *= initialSpeed;

        // Create the electron
        std::string trackFileExt{make_uuid()};
        TString trackFile{trackDir +
                          Form("/track_%s.root", trackFileExt.data())};
        ElectronTrajectoryGen traj(trackFile, field, startPos, initialVel,
                                   simStepSize, maxSimTime);
        // Check that the electron is actually trapped
        TFile fTrack(trackFile, "READ");
        TTree* tTrack = (TTree*)fTrack.Get("tree");
        double zEnd{0};
        tTrack->SetBranchAddress("zPos", &zEnd);
        tTrack->GetEntry(tTrack->GetEntries() - 1);

        if (abs(zEnd) > 10e-2) {
          cout << "Electron escaped trap. Skipping.\n";
          // Delete the track file
          gSystem->Exec("rm -f " + trackFile);
          simCounter++;
          continue;
        }
        fTrack.Close();

        // Now generate the signal
        Signal signal(trackFile, wg, lo, sampleRate);
        // Get the output voltage
        auto grV{signal.GetVITimeDomain()};

        // Delete the track file
        gSystem->Exec("rm -f " + trackFile);

        // Create dataspace for the signal
        const unsigned int DSPACE_DIM = grV->GetN();
        const unsigned int DSPACE_RANK{1};
        hsize_t dim[] = {DSPACE_DIM};
        auto dspace = new H5::DataSpace(DSPACE_RANK, dim);

        // Create dataset and write it into the file
        const std::string DATASET_NAME(GROUP_NAME + "/signal" +
                                       std::to_string(simCounter));
        auto dataset = new H5::DataSet(file->createDataSet(
            DATASET_NAME, H5::PredType::NATIVE_DOUBLE, *dspace, plist));

        // Define attributes for the dataset
        dataset
            ->createAttribute("Time step [seconds]",
                              H5::PredType::NATIVE_DOUBLE,
                              H5::DataSpace(H5S_SCALAR))
            .write(H5::PredType::NATIVE_DOUBLE, &simStepSize);
        dataset
            ->createAttribute("Rho [metres]", H5::PredType::NATIVE_DOUBLE,
                              H5::DataSpace(H5S_SCALAR))
            .write(H5::PredType::NATIVE_DOUBLE, &rho);
        dataset
            ->createAttribute("Phi [radians]", H5::PredType::NATIVE_DOUBLE,
                              H5::DataSpace(H5S_SCALAR))
            .write(H5::PredType::NATIVE_DOUBLE, &phi);
        dataset
            ->createAttribute("Pitch angle [radians]",
                              H5::PredType::NATIVE_DOUBLE,
                              H5::DataSpace(H5S_SCALAR))
            .write(H5::PredType::NATIVE_DOUBLE, &pitchAngle);
        dataset
            ->createAttribute("Kinetic energy [eV]",
                              H5::PredType::NATIVE_DOUBLE,
                              H5::DataSpace(H5S_SCALAR))
            .write(H5::PredType::NATIVE_DOUBLE, &kineticEnergy);
        dataset
            ->createAttribute("xStart [metres]", H5::PredType::NATIVE_DOUBLE,
                              H5::DataSpace(H5S_SCALAR))
            .write(H5::PredType::NATIVE_DOUBLE, &xStart);
        dataset
            ->createAttribute("yStart [metres]", H5::PredType::NATIVE_DOUBLE,
                              H5::DataSpace(H5S_SCALAR))
            .write(H5::PredType::NATIVE_DOUBLE, &yStart);
        dataset
            ->createAttribute("f_LO [Hertz]", H5::PredType::NATIVE_DOUBLE,
                              H5::DataSpace(H5S_SCALAR))
            .write(H5::PredType::NATIVE_DOUBLE, &loFreq);

        // Trap and coil config
        dataset
            ->createAttribute("r_coil [metres]", H5::PredType::NATIVE_DOUBLE,
                              H5::DataSpace(H5S_SCALAR))
            .write(H5::PredType::NATIVE_DOUBLE, &rCoil);
        dataset
            ->createAttribute("i_coil [Amps]", H5::PredType::NATIVE_DOUBLE,
                              H5::DataSpace(H5S_SCALAR))
            .write(H5::PredType::NATIVE_DOUBLE, &iCoil);
        dataset
            ->createAttribute("B_bkg [Tesla]", H5::PredType::NATIVE_DOUBLE,
                              H5::DataSpace(H5S_SCALAR))
            .write(H5::PredType::NATIVE_DOUBLE, &bkgField);
        dataset
            ->createAttribute("r_wg [metres]", H5::PredType::NATIVE_DOUBLE,
                              H5::DataSpace(H5S_SCALAR))
            .write(H5::PredType::NATIVE_DOUBLE, &wgRadius);

        // Create the buffer for writing in the voltage data
        double bufferIn[DSPACE_DIM];
        for (uint i{0}; i < grV->GetN(); i++) {
          bufferIn[i] = grV->GetPointY(i);
        }
        dataset->write(bufferIn, H5::PredType::NATIVE_DOUBLE, *dspace);

        simCounter++;
        const clock_t endEventClock = clock();
        cout << "Event generation took "
             << double(endEventClock - startEventClock) / CLOCKS_PER_SEC
             << " s\n";

        delete dataset;
        delete dspace;
      }
    }
  }

  delete group;
  delete file;

  delete field;
  delete wg;

  return 0;
}