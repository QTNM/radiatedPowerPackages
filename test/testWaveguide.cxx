/*
  testWaveguide.cxx

  Compiles to executable testWaveguide which generates validation plots for
  circular and rectangular waveguides.
*/

#include <getopt.h>
#include <unistd.h>

#include <boost/lexical_cast.hpp>
#include <boost/uuid/uuid.hpp>
#include <boost/uuid/uuid_generators.hpp>
#include <boost/uuid/uuid_io.hpp>
#include <iostream>
#include <string>

#include "BasicFunctions/BasicFunctions.h"
#include "ElectronDynamics/QTNMFields.h"
#include "ElectronDynamics/TrajectoryGen.h"
#include "H5Cpp.h"
#include "SignalProcessing/LocalOscillator.h"
#include "SignalProcessing/Signal.h"
#include "TSystem.h"
#include "Waveguides/CircularWaveguide.h"
#include "Waveguides/IWaveguide.h"
#include "Waveguides/Probe.h"
#include "Waveguides/RectangularWaveguide.h"
#include "Waveguides/WaveguideMode.h"

/// @brief Function for generating a random file string
/// @return A unique string which can be used for file names
std::string make_uuid() {
  return boost::lexical_cast<std::string>((boost::uuids::random_generator())());
}

/// @brief Function for generating data of the electric field across a
/// waveguide.
/// @param theData The data to be filled
/// @param wg The waveguide to be used
/// @param mode The mode to be used
/// @param omega The angular frequency of the wave
/// @param nPoints The number of points to sample along each axis
/// @param xMinMax The maximum value of the x axis
/// @param yMinMax The maximum value of the y axis
/// @param state The polarisation state
template <size_t rows, size_t cols>
void CalcFieldData(double (&theData)[rows][cols], rad::IWaveguide *wg,
                   rad::WaveguideMode mode, double omega, unsigned int nPoints,
                   double xMinMax = 5, double yMinMax = 5, bool state = true) {
  const double dX{2 * xMinMax / double(nPoints - 1)};
  const double dY{2 * yMinMax / double(nPoints - 1)};
  const double xMin{-xMinMax};
  const double yMin{-yMinMax};
  const double zPos{0.0};
  for (unsigned int iX{0}; iX < nPoints; ++iX) {
    double xPos{xMin + dX * double(iX)};
    for (unsigned int iY{0}; iY < nPoints; ++iY) {
      double yPos{yMin + dY * double(iY)};
      TVector3 pos(xPos, yPos, zPos);
      TVector3 eField{wg->GetModeEField(pos, mode, 1.0, omega, state)};
      theData[iX][iY] = eField.Mag();
    }
  }
}

int main(int argc, char *argv[]) {
  int opt{};
  std::string outputStr{" "};  // Location of output file
  unsigned int nPoints{40};    // Number of points to sample in the waveguide

  while ((opt = getopt(argc, argv, ":o:n:")) != -1) {
    switch (opt) {
      case 'o':
        outputStr = optarg;
        break;
      case 'n':
        nPoints = std::stoi(optarg);
        if (nPoints < 1) {
          std::cerr << "Number of points must be greater than 0.\n";
          return 1;
        }
        break;
      case ':':
        std::cerr << "Option -" << static_cast<char>(optopt)
                  << " requires an argument.\n";
        return 1;
      case '?':
        std::cerr << "Unknown option: " << static_cast<char>(optopt) << "\n";
        return 1;
    }
  }
  std::cout << "Writing output to " << outputStr << "\n";

  // Measure the collected power as a function of position
  auto fout = new H5::H5File(outputStr, H5F_ACC_TRUNC);  // Open the file
  std::string rectGroupName{"/WR42"};
  auto rectGroup = new H5::Group(fout->createGroup(rectGroupName));
  double fillValue{0};
  H5::DSetCreatPropList plist;
  plist.setFillValue(H5::PredType::NATIVE_DOUBLE, &fillValue);

  // Set up a flat background field
  const double bField{1.0};  // T
  rad::UniformField *field = new rad::UniformField(bField);

  // Define basic electron parameters
  const double eKE{18.6e3};  // eV
  const double eSpeed{rad::GetSpeedFromKE(eKE, rad::ME)};
  const TVector3 eVel(eSpeed, 0, 0);
  const double gyroradius{rad::GetGyroradius(
      eVel, field->evaluate_field_at_point(TVector3(0, 0, 0)), rad::ME)};
  const double cyclotronFreq{rad::CalcCyclotronFreq(eKE, bField)};  // Hertz
  // Calculate the free space radiated power
  const double freeSpacePower{rad::CalcLarmorPower(eKE, bField, M_PI / 2)};

  const double loadImpedance{50};  // Ohm

  // Create a waveguide with dimensions of a WR42
  const double WR42_A{10.668e-3};      // m
  const double WR42_B{4.318e-3};       // m
  const double waveguideLength{0.05};  // m
  TVector3 probePosition1(0, 0, waveguideLength / 2);
  TVector3 probePosition2(0, 0, -waveguideLength / 2);
  rad::RectangularWaveguide *wr42 =
      new rad::RectangularWaveguide(WR42_A, WR42_B, waveguideLength);

  const double xPosMinWR42{-wr42->GetLongDimension() * 0.9 / 2.0};
  const double xPosMaxWR42{wr42->GetLongDimension() * 0.9 / 2.0};

  const double simStepSize{1e-12};  // seconds
  const double simTime{1e-7};       // seconds

  const unsigned int DSPACE_DIM{nPoints};
  const int DSPACE_RANK{2};
  hsize_t dim[] = {DSPACE_DIM, 2};
  auto dspace = new H5::DataSpace(DSPACE_RANK, dim);
  const std::string DATASET_NAME(rectGroupName + "/collectedPower");
  auto dataset = new H5::DataSet(rectGroup->createDataSet(
      DATASET_NAME, H5::PredType::NATIVE_DOUBLE, *dspace, plist));

  double powerFractions[nPoints][2];
  for (unsigned int iP{0}; iP < nPoints; ++iP) {
    double xPos{xPosMinWR42 +
                (xPosMaxWR42 - xPosMinWR42) * double(iP) / double(nPoints - 1)};
    TVector3 ePos(xPos, -gyroradius, 0);

    // Generate the electron trajectory
    std::string trackFileExt{make_uuid()};
    TString trackFile{Form("track_%s.root", trackFileExt.data())};
    rad::ElectronTrajectoryGen traj(trackFile, field, ePos, eVel, simStepSize,
                                    simTime);

    // Generate the signal on the probe
    const double sampleRate{1e9};                         // Hz
    const double loFreq{cyclotronFreq - sampleRate / 4};  // Hz
    rad::LocalOscillator lo(2 * M_PI * loFreq);
    rad::Probe pr1(probePosition1, rad::WaveguideMode(1, 0, rad::kTE));
    rad::Probe pr2(probePosition2, rad::WaveguideMode(1, 0, rad::kTE));
    rad::Signal sig1(trackFile, wr42, lo, sampleRate, pr1);
    rad::Signal sig2(trackFile, wr42, lo, sampleRate, pr2);
    // Get the output voltage and calculate the average power
    auto grV1{sig1.GetVITimeDomain()};
    auto grV2{sig2.GetVITimeDomain()};
    double avgPower{0};
    for (int i{0}; i < grV1->GetN(); ++i) {
      double voltage1{grV1->GetPointY(i)};
      double voltage2{grV2->GetPointY(i)};
      avgPower += (pow(voltage1, 2) + pow(voltage2, 2)) / loadImpedance;
    }
    avgPower /= double(grV1->GetN());
    powerFractions[iP][0] = xPos;
    powerFractions[iP][1] = avgPower / freeSpacePower;

    // Delete the trajectory file
    gSystem->Exec("rm -f " + trackFile);
  }
  // Write the data to the file
  dataset->write(powerFractions, H5::PredType::NATIVE_DOUBLE, *dspace);
  delete dataset;
  delete rectGroup;

  // Now do the same for a circular waveguide
  std::string circGroupName{"/Circular"};
  auto circGroup = new H5::Group(fout->createGroup(circGroupName));
  // Define circular waveguide
  const double wgRadius = {5e-3};  // m
  auto circwg = new rad::CircularWaveguide(wgRadius, waveguideLength);
  const std::string CIRC_DATASET_NAME(circGroupName + "/collectedPower");
  auto dataset2 = new H5::DataSet(circGroup->createDataSet(
      CIRC_DATASET_NAME, H5::PredType::NATIVE_DOUBLE, *dspace, plist));

  const double xPosMinCirc{-circwg->GetInnerRadius() * 0.9};
  const double xPosMaxCirc{circwg->GetInnerRadius() * 0.9};

  for (unsigned int iP{0}; iP < nPoints; ++iP) {
    double xPos{xPosMinCirc +
                (xPosMaxCirc - xPosMinCirc) * double(iP) / double(nPoints - 1)};
    TVector3 ePos(xPos, -gyroradius, 0);

    // Generate the electron trajectory
    std::string trackFileExt{make_uuid()};
    TString trackFile{Form("track_%s.root", trackFileExt.data())};
    rad::ElectronTrajectoryGen traj(trackFile, field, ePos, eVel, simStepSize,
                                    simTime);

    // Generate the signal on the probe
    const double sampleRate{1e9};                         // Hz
    const double loFreq{cyclotronFreq - sampleRate / 4};  // Hz
    rad::LocalOscillator lo(2 * M_PI * loFreq);

    rad::Probe pr1_1(probePosition1, rad::WaveguideMode(1, 1, rad::kTE), true);
    rad::Probe pr1_2(probePosition1, rad::WaveguideMode(1, 1, rad::kTE), false);
    rad::Probe pr2_1(probePosition2, rad::WaveguideMode(1, 1, rad::kTE), true);
    rad::Probe pr2_2(probePosition2, rad::WaveguideMode(1, 1, rad::kTE), false);

    rad::Signal sig1_1(trackFile, circwg, lo, sampleRate, pr1_1);
    rad::Signal sig1_2(trackFile, circwg, lo, sampleRate, pr1_2);
    rad::Signal sig2_1(trackFile, circwg, lo, sampleRate, pr2_1);
    rad::Signal sig2_2(trackFile, circwg, lo, sampleRate, pr2_2);
    // Get the output voltage and calculate the average power
    auto grV1_1{sig1_1.GetVITimeDomain()};
    auto grV1_2{sig1_2.GetVITimeDomain()};
    auto grV2_1{sig2_1.GetVITimeDomain()};
    auto grV2_2{sig2_2.GetVITimeDomain()};
    double avgPower{0};
    for (int i{0}; i < grV1_1->GetN(); ++i) {
      double voltage1_1{grV1_1->GetPointY(i)};
      double voltage1_2{grV1_2->GetPointY(i)};
      double voltage2_1{grV2_1->GetPointY(i)};
      double voltage2_2{grV2_2->GetPointY(i)};
      avgPower += (pow(voltage1_1, 2) + pow(voltage1_2, 2) +
                   pow(voltage2_1, 2) + pow(voltage2_2, 2)) /
                  loadImpedance;
    }
    avgPower /= double(grV1_1->GetN());
    powerFractions[iP][0] = xPos;
    powerFractions[iP][1] = avgPower / freeSpacePower;

    // Delete the trajectory file
    gSystem->Exec("rm -f " + trackFile);
  }

  // Write the data to the file
  dataset2->write(powerFractions, H5::PredType::NATIVE_DOUBLE, *dspace);
  delete dataset2;

  delete dspace;
  delete circGroup;

  /////////////////////////////////////////////////////////////
  /////////////// Field pattern generation ////////////////////
  /////////////////////////////////////////////////////////////
  // Create a new group for the field data
  std::string fieldGroupName{"/FieldData"};
  auto fieldGroup = new H5::Group(fout->createGroup(fieldGroupName));
  // Set up the data space
  const unsigned int nFieldPoints{41};
  const double xMinMax{5e-3};
  const double yMinMax{5e-3};
  double circWgFieldMags_P[nFieldPoints][nFieldPoints];
  double circWgFieldMags_M[nFieldPoints][nFieldPoints];
  double rectWgFieldMags[nFieldPoints][nFieldPoints];
  // Set up the data space
  hsize_t dimField[] = {nFieldPoints, nFieldPoints};
  auto dspaceField = new H5::DataSpace(2, dimField);
  const std::string CIRC_FIELD_P_DATASET_NAME(fieldGroupName +
                                              "/circFieldDataTE11_P");
  auto datasetCircP = new H5::DataSet(fieldGroup->createDataSet(
      CIRC_FIELD_P_DATASET_NAME, H5::PredType::NATIVE_DOUBLE, *dspaceField,
      plist));
  const std::string CIRC_FIELD_M_DATASET_NAME(fieldGroupName +
                                              "/circFieldDataTE11_M");
  auto datasetCircM = new H5::DataSet(fieldGroup->createDataSet(
      CIRC_FIELD_M_DATASET_NAME, H5::PredType::NATIVE_DOUBLE, *dspaceField,
      plist));

  const std::string RECT_FIELD_DATASET_NAME(fieldGroupName +
                                            "/rectFieldDataTE10");
  auto datasetRect = new H5::DataSet(fieldGroup->createDataSet(
      RECT_FIELD_DATASET_NAME, H5::PredType::NATIVE_DOUBLE, *dspaceField,
      plist));

  // Get the modes to be tested
  rad::WaveguideMode circMode(1, 1, rad::kTE);
  rad::WaveguideMode rectMode(1, 0, rad::kTE);
  const double omega{2 * M_PI * cyclotronFreq};
  // Get the field data
  CalcFieldData(circWgFieldMags_P, circwg, circMode, omega, nFieldPoints,
                xMinMax, yMinMax, true);
  CalcFieldData(circWgFieldMags_M, circwg, circMode, omega, nFieldPoints,
                xMinMax, yMinMax, false);
  CalcFieldData(rectWgFieldMags, wr42, rectMode, omega, nFieldPoints, xMinMax,
                yMinMax);
  // Write the data to the file
  datasetCircP->write(circWgFieldMags_P, H5::PredType::NATIVE_DOUBLE,
                      *dspaceField);
  datasetCircM->write(circWgFieldMags_M, H5::PredType::NATIVE_DOUBLE,
                      *dspaceField);
  datasetRect->write(rectWgFieldMags, H5::PredType::NATIVE_DOUBLE,
                     *dspaceField);
  delete datasetCircP;
  delete datasetCircM;
  delete datasetRect;
  delete fieldGroup;
  delete dspaceField;

  delete wr42;
  delete circwg;

  delete fout;
  return 0;
}