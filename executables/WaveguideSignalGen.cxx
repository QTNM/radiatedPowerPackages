/*
  WaveguideSignalGen.cxx

  Use this to generate ensembles of signals in waveguides
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
#include "SignalProcessing/LocalOscillator.h"
#include "SignalProcessing/Signal.h"
#include "TFile.h"
#include "TGraph.h"
#include "TSystem.h"
#include "TTree.h"
#include "TVector3.h"
#include "Waveguides/CircularWaveguide.h"

using namespace rad;
using std::cout;
using std::endl;

bool GenerateElectron(TString file, TVector3 pos, TVector3 vel,
                      BaseField *field, long double stepSize, double simTime) {
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

  unsigned long nSteps{1};
  while (time <= simTime && isTrapped) {
    time = (long double)nSteps * stepSize;
    if (std::fmod(time, 5e-6) < stepSize) {
      std::cout << "Simulated " << time * 1e6 << " us\n";
    }

    // Advance the electron
    auto outputStep = solver.advance_step(stepSize, p, v);
    p = std::get<0>(outputStep);
    if (abs(p.Z()) > 15e-2) {
      std::cout << "Electron escaped after " << time * 1e6 << " us\n";
      isTrapped = false;
      break;
    }

    v = std::get<1>(outputStep);
    acc = solver.acc(p, v);

    // Write to tree
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
    nSteps++;
  }

  fT.cd();
  tree.Write("", TObject::kOverwrite);
  fT.Close();
  return isTrapped;
}

/// @brief Function for generating a random file string
/// @return A unique string which can be used for file names
std::string make_uuid() {
  return boost::lexical_cast<std::string>((boost::uuids::random_generator())());
}

/// @brief Main executable
int main(int argc, char *argv[]) {
  // Start a clock
  const clock_t startTotalClock{clock()};

  int opt{};
  double maxSimTime{1e-3};         // seconds
  std::string outputStemStr{" "};  // Directory in which to store files
  unsigned int nEvents{1000};      // Number of electrons to attempt to generate
  double maxRunTime{-1};           // Maximum run time in seconds

  while ((opt = getopt(argc, argv, ":o:t:n:r:")) != -1) {
    switch (opt) {
      case 'o':
        outputStemStr = optarg;
        break;

      case 't':
        maxSimTime = std::stod(optarg);
        break;

      case 'n':
        nEvents = std::stoi(optarg);
        break;

      case 'r':
        maxRunTime = std::stod(optarg);
        break;

      case ':':
        std::cout << "Option needs a value\n";
        break;

      case '?':
        std::cout << "Unknown option: " << optopt << std::endl;
        break;
    }
  }

  // Do some output checking
  if (maxSimTime <= 0) {
    cout << "Selected an invalid electron propagation time of " << maxSimTime
         << " seconds. Exiting.\n";
    exit(1);
  }

  cout << "Attempting to write out to directory: " << outputStemStr << endl;
  cout << "Simulating electron for " << maxSimTime * 1e6 << " us\n";
  cout << "Attempting to generate " << nEvents << " electrons\n";

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

  // Convert to TString
  TString outputStem{outputStemStr};

  // Create output file for end products
  std::string outfileExt{make_uuid()};
  TString outfile{outputStem + Form("/out_%s.root", outfileExt.data())};
  TFile fout(outfile, "recreate");

  // Also want to try writing this to a HDF5 file
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

  // Constant electron kinematics
  const double eKE{18.6e3};
  const double eSpeed{GetSpeedFromKE(eKE, ME)};  // m/s

  // Define the field. Do a harmonic trap.
  const double rCoil{15e-3};                        // metres
  const double trapDepth{3.4e-3};                   // Tesla
  const double iCoil{2 * trapDepth * rCoil / MU0};  // Amps
  // Magnitude of background field
  const double bkgField{0.7};  // Tesla
  auto field = new HarmonicField(rCoil, iCoil, bkgField);
  TVector3 centralField{field->evaluate_field_at_point(TVector3(0, 0, 0))};
  const double centralFieldMag{centralField.Mag()};
  const double centralCycFreq{CalcCyclotronFreq(eKE, centralFieldMag)};

  // Do about 15 time steps per cyclotron orbit
  const double simStepSize{1 / (15 * centralCycFreq)};  // seconds

  // Now define the waveguide that we want to collect our signal with
  const double wgLength{10e-2};           // metres
  const double wgRadius{7.14e-3};         // metres
  TVector3 probePos{0, 0, wgLength / 2};  // Place probe at end of guide
  auto wg = new CircularWaveguide(wgRadius, wgLength, probePos);

  // Random number stuff
  std::random_device rd;
  std::mt19937 gen(rd());

  for (uint iEv{0}; iEv < nEvents; iEv++) {
    const double rGenMax{6e-3};  // metres
    const double zGen{0};        // metres
    const clock_t startEventClock{clock()};

    // Create the uniform distribution
    std::uniform_real_distribution<double> uni1(0, 1);
    // Generate the initial electron position
    const double rGen{rGenMax * sqrt(uni1(gen))};
    const double thetaPosGen{uni1(gen) * 2 * M_PI};
    TVector3 initPos(rGen * cos(thetaPosGen), rGen * sin(thetaPosGen), zGen);
    double xStart{initPos.X()};
    double yStart{initPos.Y()};
    double zStart{initPos.Z()};

    // Generate the initial velocity
    const double phiVelGen{uni1(gen) * 2 * M_PI};
    const double thetaVelGen{uni1(gen) * M_PI};
    TVector3 initVel(eSpeed * cos(phiVelGen) * sin(thetaVelGen),
                     eSpeed * sin(phiVelGen) * sin(thetaVelGen),
                     eSpeed * cos(thetaVelGen));

    // Calculate pitch angle
    const double pitchAngleRad{abs(atan(initVel.Perp() / initVel.Z()))};
    const double pitchAngleDeg{pitchAngleRad * 180 / M_PI};
    cout << "Event " << iEv << ": Pitch angle = " << pitchAngleDeg
         << " degrees\n";

    // Actually generate the trajectory
    std::string trackFileExt{make_uuid()};
    TString trackFile{outputStem + Form("/track_%s.root", trackFileExt.data())};
    bool isTrapped = GenerateElectron(trackFile, initPos, initVel, field,
                                      simStepSize, maxSimTime);

    // Only do signal processing if our electron is trapped
    if (isTrapped) {
      // Signal processing stuff
      const double sampleRate{1e9};  // Samples per second
      const double loFreq{centralCycFreq - sampleRate / 4};  // Hz
      // Define local oscillator which should put us in the middle
      // of the sample space
      LocalOscillator lo(2 * M_PI * loFreq);
      // Do noiseless for now since it's not normalised properly
      Signal sig(trackFile, wg, lo, sampleRate);

      // We want to have the true carrier frequency as an attribute
      // so reopen the track file
      TFile fTrack(trackFile, "read");
      auto tr = (TTree *)fTrack.Get("tree");
      double time{};
      double xPos{}, yPos{}, zPos{};
      tr->SetBranchAddress("time", &time);
      tr->SetBranchAddress("xPos", &xPos);
      tr->SetBranchAddress("yPos", &yPos);
      tr->SetBranchAddress("zPos", &zPos);
      double zCrossTimes[3];
      uint zCrossings{0};
      tr->GetEntry(0);
      double oldZ{zPos};
      for (int iE{1}; iE < tr->GetEntries(); iE++) {
        tr->GetEntry(iE);

        if ((oldZ > 0 && zPos < 0) || (oldZ < 0 && zPos > 0)) {
          zCrossTimes[zCrossings] = time;
          zCrossings++;
          if (zCrossings > 2) break;
        }
        oldZ = zPos;
      }

      double bMean{0};
      double zMax{-DBL_MAX};
      uint nFieldPnts{0};
      for (int iE{0}; iE < tr->GetEntries(); iE++) {
        tr->GetEntry(iE);
        TVector3 ePos(xPos, yPos, zPos);
        if (time > 1e-6) break;

        bMean += field->evaluate_field_magnitude(ePos);
        nFieldPnts++;

	if (zPos > zMax) 
	  zMax = zPos;
      }
      bMean /= double(nFieldPnts);

      double startFreq{CalcCyclotronFreq(eKE, bMean)};
      double startFreqDM{CalcCyclotronFreq(eKE, bMean) - loFreq};
      double axFreq{1 / (zCrossTimes[2] - zCrossTimes[0])};
      cout << "Initial frequency = " << startFreq / 1e6
           << " MHz\tAxial frequency = " << axFreq / 1e6 << " MHz\n";
      delete tr;
      fTrack.Close();

      // Write the TGraph actually containing the information about the signal
      // to the output file
      auto grV{sig.GetVITimeDomain()};
      fout.cd();
      grV->Write(Form("grV%d", iEv));
      auto grVSpec{MakePowerSpectrumPeriodogram(grV)};
      ScaleGraph(grVSpec, 1e15);
      grVSpec->SetTitle(
          Form("%.1f MHz, #theta = %.1f degrees, r_{i} = %.1f mm; Frequency "
               "[Hz]; Power [fW]",
               startFreq / 1e6, pitchAngleDeg, rGen * 1e3));
      grVSpec->Write(Form("grVSpec%d", iEv));

      // Create dataspace for dataset
      const unsigned int DSPACE_DIM = grV->GetN();
      cout << "Time series is " << DSPACE_DIM << " entries long\n";
      const unsigned int DSPACE_RANK{1};
      hsize_t dim[] = {DSPACE_DIM};
      auto dspace = new H5::DataSpace(DSPACE_RANK, dim);

      // Create dataset and write it into the file
      const std::string DATASET_NAME(GROUP_NAME + "/signal" +
                                     std::to_string(iEv));
      auto dataset = new H5::DataSet(file->createDataSet(
          DATASET_NAME, H5::PredType::NATIVE_DOUBLE, *dspace, plist));

      // Define attributes
      // First create the attribute for the time step
      H5::DataSpace timeStepSpace(H5S_SCALAR);
      H5::Attribute timeStepAttr{dataset->createAttribute(
          "Time step (seconds)", H5::PredType::NATIVE_DOUBLE, timeStepSpace)};
      double timeStep{grV->GetPointX(1) - grV->GetPointX(0)};
      timeStepAttr.write(H5::PredType::NATIVE_DOUBLE, &timeStep);
      // Now do the pitch angle
      H5::DataSpace pitchSpace(H5S_SCALAR);
      H5::Attribute pitchAttr{dataset->createAttribute(
          "Pitch angle (degrees)", H5::PredType::NATIVE_DOUBLE, pitchSpace)};
      pitchAttr.write(H5::PredType::NATIVE_DOUBLE, &pitchAngleDeg);
      // And the radial offset
      H5::DataSpace radSpace(H5S_SCALAR);
      H5::Attribute radAttr{dataset->createAttribute(
          "Radial offset (metres)", H5::PredType::NATIVE_DOUBLE, radSpace)};
      radAttr.write(H5::PredType::NATIVE_DOUBLE, &rGen);

      // Start x
      H5::DataSpace xStartSpc(H5S_SCALAR);
      H5::Attribute xStartAttr{dataset->createAttribute(
          "x start (metres)", H5::PredType::NATIVE_DOUBLE, xStartSpc)};
      xStartAttr.write(H5::PredType::NATIVE_DOUBLE, &xStart);
      // Start y
      H5::DataSpace yStartSpc(H5S_SCALAR);
      H5::Attribute yStartAttr{dataset->createAttribute(
          "y start (metres)", H5::PredType::NATIVE_DOUBLE, yStartSpc)};
      yStartAttr.write(H5::PredType::NATIVE_DOUBLE, &yStart);
      // Start z
      H5::DataSpace zStartSpc(H5S_SCALAR);
      H5::Attribute zStartAttr{dataset->createAttribute(
          "z start (metres)", H5::PredType::NATIVE_DOUBLE, zStartSpc)};
      zStartAttr.write(H5::PredType::NATIVE_DOUBLE, &zStart);

      // Start frequency
      H5::DataSpace startFreqSpc(H5S_SCALAR);
      H5::Attribute startFreqAttr{
          dataset->createAttribute("Start frequency (Hertz)",
                                   H5::PredType::NATIVE_DOUBLE, startFreqSpc)};
      startFreqAttr.write(H5::PredType::NATIVE_DOUBLE, &startFreq);
      // Start frequency (downmixed)
      H5::DataSpace startFreqDMSpc(H5S_SCALAR);
      H5::Attribute startFreqDMAttr{
          dataset->createAttribute("Start frequency, downmixed (Hertz)",
                                   H5::PredType::NATIVE_DOUBLE, startFreqDMSpc)};
      startFreqDMAttr.write(H5::PredType::NATIVE_DOUBLE, &startFreqDM);
      // Axial frequency
      H5::DataSpace axFreqSpc(H5S_SCALAR);
      H5::Attribute axFreqAttr{dataset->createAttribute(
          "Axial frequency (Hertz)", H5::PredType::NATIVE_DOUBLE, axFreqSpc)};
      axFreqAttr.write(H5::PredType::NATIVE_DOUBLE, &axFreq);

      // Max displacement
      H5::DataSpace zMaxSpc(H5S_SCALAR);
      H5::Attribute zMaxAttr{
	dataset->createAttribute("z_max (metres)", H5::PredType::NATIVE_DOUBLE, 
				 zMaxSpc)};
      zMaxAttr.write(H5::PredType::NATIVE_DOUBLE, &zMax);


      // Create the buffer for writing in data
      double bufferIn[DSPACE_DIM];
      for (uint i{0}; i < grV->GetN(); i++) {
        bufferIn[i] = grV->GetPointY(i);
      }
      dataset->write(bufferIn, H5::PredType::NATIVE_DOUBLE, *dspace);

      // Close the dataset
      delete dataset;
      // Close the dataspace
      delete dspace;

      delete grVSpec;
    }

    // The track files can get pretty large so it's best to delete them after
    // you're done using them
    gSystem->Exec("rm -f " + trackFile);

    const clock_t endEventClock{clock()};
    const double eventTime{double(endEventClock - startEventClock) /
                           CLOCKS_PER_SEC};
    const double runTime{double(endEventClock - startTotalClock) /
                         CLOCKS_PER_SEC};
    cout << "Event took " << eventTime
         << " seconds to generate\tTotal time = " << runTime << "s\n\n";

    // Stop generating events if we're getting close to the time limit
    // Assume we're going to be taking roughly 4200s to generate each event
    if (maxRunTime > 0 & runTime > maxRunTime - 4400) break;
  }

  fout.Close();

  // Close the group
  delete group;
  // Close the file
  delete file;

  // Clock at the end of the executable
  const clock_t endTotalClock{clock()};
  double eventTime{double(endTotalClock - startTotalClock) / CLOCKS_PER_SEC};
  cout << "Took " << eventTime
       << " seconds to generate and process all events\n";

  delete field;
  delete wg;
  return 0;
}
