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
  TVector3 startPos;  // metres
  TVector3 startVel;  // m/s

  // Default constructor
  ElectronInfo() {
    startTime = 0;
    startKE = -1;
    startF = -1;
    axialF = -1;
    pitchAngle = DBL_MAX;
    zMax = -DBL_MAX;
    startPos = TVector3(0, 0, 0);
    startVel = TVector3(0, 0, 0);
  }

  void Reset() {
    startTime = 0;
    startKE = -1;
    startF = -1;
    axialF = -1;
    pitchAngle = DBL_MAX;
    zMax = -DBL_MAX;
    startPos = TVector3(0, 0, 0);
    startVel = TVector3(0, 0, 0);
  }
};

double GenerateElectron(TString file, TVector3 pos, TVector3 vel,
                        BaseField *field, double gasDensity,
                        long double stepSize, double maxSimTime,
                        std::vector<ElectronInfo> &ei, BaseField *bField) {
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

  // Now calculate when the first scatter will happen
  double gamma{1 / sqrt(1 - pow(vel.Mag() / TMath::C(), 2))};
  double ke{(gamma - 1) * ME_EV};
  ElasticScatter scatterElastic(ke);
  InelasticScatter scatterInelastic(ke);
  double elasticXSec{scatterElastic.GetTotalXSec()};
  double inelasticXSec{scatterInelastic.GetTotalXSec()};
  double totalXSec{elasticXSec + inelasticXSec};
  // Calculate the mean free path
  double lambdaStep{1 / (gasDensity * totalXSec)};
  std::exponential_distribution<double> pathDist(1 / lambdaStep);
  double pathLength{pathDist(gen)};
  double pathTimeStep{pathLength / vel.Mag()};
  double nextScatterTime{pathTimeStep};
  cout << "Scattering after " << pathTimeStep * 1e6 << " us\n";

  // Set some information about the electron creation
  ElectronInfo eiTemp;
  eiTemp.startTime = time;
  eiTemp.startKE = ke;
  eiTemp.startPos = pos;
  eiTemp.startVel = vel;

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
    time = (long double)nSteps * stepSize;
    currentInfoWriteTime += stepSize;
    if (std::fmod(time, 5e-6) < stepSize) {
      std::cout << "Simulated " << time * 1e6 << " us\n";
    }

    // Advance the electron
    auto outputStep = solver.advance_step(stepSize, p, v);
    p = std::get<0>(outputStep);
    v = std::get<1>(outputStep);
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

    if (time < nextScatterTime) {
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
        gamma = 1 / sqrt(1 - pow(v.Mag() / TMath::C(), 2));
        ke = (gamma - 1) * ME_EV;
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
    } else {
      // Rest some variables used for calculating truth info
      currentInfoWriteTime = 0;
      zCrossings = 0;
      oldZ = zPos;
      bMean = 0;
      nFieldPoints = 0;
      eiTemp.Reset();
      haveSavedInfo = false;

      // Time to do a scattering calculation
      gamma = 1 / sqrt(1 - pow(v.Mag() / TMath::C(), 2));
      ke = (gamma - 1) * ME_EV;

      eiTemp.startTime = time;
      eiTemp.startKE = ke;
      eiTemp.startPos = p;

      // Recalculate cross-sections based on the energies
      scatterElastic.SetIncidentKE(ke);
      scatterInelastic.SetIncidentKE(ke);
      elasticXSec = scatterElastic.GetTotalXSec();
      inelasticXSec = scatterInelastic.GetTotalXSec();
      totalXSec = elasticXSec + inelasticXSec;
      // Figure out if this is an elastic or inelastic scatter
      double scatterAngle{0};
      std::uniform_real_distribution<double> uni(0, 1);
      if (uni(gen) < elasticXSec / totalXSec) {
        // Elastic scatter
        cout << "Elastic scatter\n";
        scatterAngle = scatterElastic.GetRandomScatteringAngle();
        v = scatterElastic.GetScatteredVector(v, ke, scatterAngle);
      } else {
        // Inelastic scatter
        cout << "Inelastic scatter\n";
        // Get the energy of the electron after the scatter
        double W{scatterInelastic.GetRandomW()};
        scatterAngle = scatterInelastic.GetRandomTheta(W);
        // Now calculate the scattered vector
        ke = W;
        v = scatterInelastic.GetScatteredVector(v, ke, scatterAngle);
      }
      double pitchAngleRad{abs(atan(v.Perp() / v.Z()))};
      double pitchAngleDeg{pitchAngleRad * 180 / M_PI};
      cout << "Scattering angle = " << scatterAngle * 180 / M_PI
           << " degrees\tNew KE = " << ke
           << " keV\tNew pitch angle = " << pitchAngleDeg << " degrees\n";
      eiTemp.startVel = v;

      // Calculate the new acceleration
      acc = solver.acc(p, v);
      // Fill the tree
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

      // Get the next scatter time
      scatterElastic.SetIncidentKE(ke);
      scatterInelastic.SetIncidentKE(ke);
      elasticXSec = scatterElastic.GetTotalXSec();
      inelasticXSec = scatterInelastic.GetTotalXSec();
      totalXSec = elasticXSec + inelasticXSec;
      lambdaStep = 1 / (gasDensity * totalXSec);
      std::exponential_distribution<double> pathDistNext(1 / lambdaStep);
      pathLength = pathDistNext(gen);
      pathTimeStep = pathLength / v.Mag();
      nextScatterTime += pathTimeStep;
      cout << "Next scatter at " << nextScatterTime * 1e6 << " us\n";
    }
    nSteps++;
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
  double tritiumGasDensity{1e18};  // m^-3
  double maxSimTime{10e-3};        // seconds

  while ((opt = getopt(argc, argv, ":o:n:d:t:")) != -1) {
    switch (opt) {
      case 'o':
        outputStemStr = optarg;
        break;
      case 'n':
        nEvents = boost::lexical_cast<unsigned int>(optarg);
        break;
      case 'd':
        tritiumGasDensity = boost::lexical_cast<double>(optarg);
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
    }
  }

  cout << "Writing output files to " << outputStemStr << "\n";
  cout << "Attempting to generate " << nEvents << " events.\n";
  cout << "Maximum simulation time = " << maxSimTime * 1e6 << " us\n";
  cout << "Tritium gas density = " << tritiumGasDensity << " m^-3\n";

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

  /*
  TGraph *grField = new TGraph();
  setGraphAttr(grField);
  for (int i{0}; i < 400; i++) {
    double z{-0.1 + double(i) * 0.2 / 399.0};
    grField->SetPoint(i, z, field->evaluate_field_magnitude(TVector3(0, 0, z)));
  }
  rootFile.cd();
  grField->Write("Field");
  */

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

  // Define some simulation parameters
  uint nGenerated{0};
  const double rhoGenMax{5e-3};     // metres
  const double energyWindow{1000};  // eV
  const double zLimit{wgLength / 8};
  while (nGenerated < nEvents) {
    const clock_t startEventClock{clock()};

    // Create a uniform distribution
    std::uniform_real_distribution<double> uni1(0, 1);
    // Generate the initial electron position
    double zGen{-zLimit + uni1(gen) * 2 * zLimit};
    const double rhoGen{rhoGenMax * sqrt(uni1(gen))};
    const double phiPosGen{uni1(gen) * 2 * M_PI};
    TVector3 initPos(rhoGen * cos(phiPosGen), rhoGen * sin(phiPosGen), zGen);
    double xStart{initPos.X()};
    double yStart{initPos.Y()};
    double zStart{initPos.Z()};

    // Generate the initial direction
    const double phiVelGen{uni1(gen) * 2 * M_PI};
    const double thetaVelGen{uni1(gen) * M_PI};

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
    cout << "Event " << nGenerated << ": Energy = " << initialKE
         << " eV\tPitch angle = " << pitchAngleDeg << " degrees\n";

    // Now generate the trajectory
    std::string trackFileExt{make_uuid()};
    TString trackFile{outputDir + Form("/track_%s.root", trackFileExt.data())};

    InelasticScatter scatterInel(initialKE);
    cout << "Inelastic cross-section = " << scatterInel.GetTotalXSec()
         << " m^2\n";

    std::vector<ElectronInfo> eiVec{};
    double simTime{GenerateElectron(trackFile, initPos, initVel, field,
                                    tritiumGasDensity, simStepSize, maxSimTime,
                                    eiVec, field)};
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

    if (simTime > 1e-6) {
      const clock_t endPropClock{clock()};
      const double propTime{(endPropClock - startEventClock) /
                            (double)CLOCKS_PER_SEC};
      cout << "Event " << nGenerated << " took " << propTime
           << " seconds to propagate electron.\n";

      // Now generate the signal
      Signal signal(trackFile, wg, lo, sampleRate);
      // Get the output voltage
      auto grV{signal.GetVITimeDomain()};

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

      // Create dataspace for dataset
      const unsigned int DSPACE_DIM = grV->GetN();
      cout << "Time series is " << DSPACE_DIM << " entries long\n";
      const unsigned int DSPACE_RANK{1};
      hsize_t dim[] = {DSPACE_DIM};
      auto dspace = new H5::DataSpace(DSPACE_RANK, dim);

      // Create dataset and write it into the file
      const std::string DATASET_NAME(GROUP_NAME + "/signal");
      auto dataset = new H5::DataSet(file->createDataSet(
          DATASET_NAME, H5::PredType::NATIVE_DOUBLE, *dspace, plist));
      // Define attributes
      // First create the attribute for the time step
      H5::DataSpace timeStepSpace(H5S_SCALAR);
      H5::Attribute timeStepAttr{dataset->createAttribute(
          "Time step [seconds]", H5::PredType::NATIVE_DOUBLE, timeStepSpace)};
      double timeStep{grV->GetPointX(1) - grV->GetPointX(0)};
      timeStepAttr.write(H5::PredType::NATIVE_DOUBLE, &timeStep);
      // Local oscillator frequency
      H5::DataSpace loFreqSpc(H5S_SCALAR);
      H5::Attribute loFreqAttr{dataset->createAttribute(
          "LO frequency [Hertz]", H5::PredType::NATIVE_DOUBLE, loFreqSpc)};
      loFreqAttr.write(H5::PredType::NATIVE_DOUBLE, &loFreq);
      // Write the metadata for the trap config
      H5::DataSpace rCoilSpc(H5S_SCALAR);
      H5::Attribute rCoilAttr{dataset->createAttribute(
          "r_coil [metres]", H5::PredType::NATIVE_DOUBLE, rCoilSpc)};
      rCoilAttr.write(H5::PredType::NATIVE_DOUBLE, &rCoil);
      H5::DataSpace iCoilSpc(H5S_SCALAR);
      H5::Attribute iCoilAttr{dataset->createAttribute(
          "i_coil [Amps]", H5::PredType::NATIVE_DOUBLE, iCoilSpc)};
      iCoilAttr.write(H5::PredType::NATIVE_DOUBLE, &iCoil);
      H5::DataSpace bkgFieldSpc(H5S_SCALAR);
      H5::Attribute bkgFieldAttr{dataset->createAttribute(
          "B_bkg [Tesla]", H5::PredType::NATIVE_DOUBLE, bkgFieldSpc)};
      bkgFieldAttr.write(H5::PredType::NATIVE_DOUBLE, &bkgField);
      // Write the metadata for the waveguide dimenstions as well
      H5::DataSpace rWgSpc(H5S_SCALAR);
      H5::Attribute rWgAttr{dataset->createAttribute(
          "r_wg [metres]", H5::PredType::NATIVE_DOUBLE, rWgSpc)};
      rWgAttr.write(H5::PredType::NATIVE_DOUBLE, &wgRadius);
      // Metadata for gas density
      H5::DataSpace gasDensitySpc(H5S_SCALAR);
      H5::Attribute gasDensityAttr{dataset->createAttribute(
          "Gas density [m^-3]", H5::PredType::NATIVE_DOUBLE, gasDensitySpc)};
      gasDensityAttr.write(H5::PredType::NATIVE_DOUBLE, &tritiumGasDensity);

      // Now create the attributes for each scattering event
      const unsigned int SCATTERSPACE_DIM = eiVec.size();
      hsize_t dimScatter[] = {SCATTERSPACE_DIM};
      H5::DataSpace scatterSpace(1, dimScatter);
      H5::Attribute timeAttr{dataset->createAttribute(
          "Time [seconds]", H5::PredType::NATIVE_DOUBLE, scatterSpace)};
      double timeBuffer[SCATTERSPACE_DIM];
      H5::Attribute startEAttr{dataset->createAttribute(
          "Start E [eV]", H5::PredType::NATIVE_DOUBLE, scatterSpace)};
      double startEBuffer[SCATTERSPACE_DIM];
      H5::Attribute startFAttr{dataset->createAttribute(
          "Start F [Hertz]", H5::PredType::NATIVE_DOUBLE, scatterSpace)};
      double startFBuffer[SCATTERSPACE_DIM];
      H5::Attribute startFDMAttr{
          dataset->createAttribute("Start F, downmixed [Hertz]",
                                   H5::PredType::NATIVE_DOUBLE, scatterSpace)};
      double startFDMBuffer[SCATTERSPACE_DIM];
      H5::Attribute axialFAttr{dataset->createAttribute(
          "Axial F [Hertz]", H5::PredType::NATIVE_DOUBLE, scatterSpace)};
      double axialFBuffer[SCATTERSPACE_DIM];
      H5::Attribute zMaxAttr{dataset->createAttribute(
          "z_max [metres]", H5::PredType::NATIVE_DOUBLE, scatterSpace)};
      double zMaxBuffer[SCATTERSPACE_DIM];
      H5::Attribute pitchAngleAttr{dataset->createAttribute(
          "Pitch angle [degrees]", H5::PredType::NATIVE_DOUBLE, scatterSpace)};
      double pitchAngleBuffer[SCATTERSPACE_DIM];

      H5::Attribute xpStartAttr{dataset->createAttribute(
          "x_start [metres]", H5::PredType::NATIVE_DOUBLE, scatterSpace)};
      double xpStartBuffer[SCATTERSPACE_DIM];
      H5::Attribute ypStartAttr{dataset->createAttribute(
          "y_start [metres]", H5::PredType::NATIVE_DOUBLE, scatterSpace)};
      double ypStartBuffer[SCATTERSPACE_DIM];
      H5::Attribute zpStartAttr{dataset->createAttribute(
          "z_start [metres]", H5::PredType::NATIVE_DOUBLE, scatterSpace)};
      double zpStartBuffer[SCATTERSPACE_DIM];

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
      }
      timeAttr.write(H5::PredType::NATIVE_DOUBLE, &timeBuffer);
      startEAttr.write(H5::PredType::NATIVE_DOUBLE, &startEBuffer);
      startFAttr.write(H5::PredType::NATIVE_DOUBLE, &startFBuffer);
      startFDMAttr.write(H5::PredType::NATIVE_DOUBLE, &startFDMBuffer);
      axialFAttr.write(H5::PredType::NATIVE_DOUBLE, &axialFBuffer);
      zMaxAttr.write(H5::PredType::NATIVE_DOUBLE, &zMaxBuffer);
      pitchAngleAttr.write(H5::PredType::NATIVE_DOUBLE, &pitchAngleBuffer);
      xpStartAttr.write(H5::PredType::NATIVE_DOUBLE, &xpStartBuffer);
      ypStartAttr.write(H5::PredType::NATIVE_DOUBLE, &ypStartBuffer);
      zpStartAttr.write(H5::PredType::NATIVE_DOUBLE, &zpStartBuffer);

      // Create the buffer for writing in the voltage data
      double bufferIn[DSPACE_DIM];
      for (uint i{0}; i < grV->GetN(); i++) {
        bufferIn[i] = grV->GetPointY(i);
      }
      dataset->write(bufferIn, H5::PredType::NATIVE_DOUBLE, *dspace);

      // Close the dataset and dataspace
      delete dataset;
      delete dspace;
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