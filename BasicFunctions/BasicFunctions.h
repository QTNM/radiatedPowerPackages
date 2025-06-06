// BasicFunctions.h
#ifndef BASIC_FUNCTIONS_H
#define BASIC_FUNCTIONS_H

#include <vector>

#include "BasicFunctions/Constants.h"
#include "BasicFunctions/FFTWComplex.h"
#include "BasicFunctions/FourierTransforms.h"
#include "Math/Point3D.h"
#include "Math/Vector3D.h"
#include "TGraph.h"
#include "TH1.h"
#include "TH2.h"
#include "TVector3.h"

namespace rad {
void setGraphAttr(TGraph *gr);

void setGraphAttr(std::unique_ptr<TGraph> &gr);

/// @brief Formatter for TGraph references
/// @param gr Graph to be formatted
void setGraphAttr(TGraph &gr);

void SetHistAttr(TH1 *h);

void SetHistAttr(TH1 &h);

void SetHistAttr(std::unique_ptr<TH1D> &h);

void SetHistAttr(TH2 *h);

void SetHistAttr(std::unique_ptr<TH2> &h);

double CalcAeHertzianDipole(const double wavelength,
                            const ROOT::Math::XYZVector dipoleDir,
                            const ROOT::Math::XYZPoint ePosition,
                            const ROOT::Math::XYZPoint antennaPoint);

double CalcAlHertzianDipole(const double wavelength,
                            const ROOT::Math::XYZVector dipoleDir,
                            const ROOT::Math::XYZPoint ePosition,
                            const ROOT::Math::XYZPoint antennaPoint);

double CalcRetardedTime(const ROOT::Math::XYZPoint fieldPoint,
                        const ROOT::Math::XYZPoint ePosition,
                        const double labTime);

/// For a given retarded time, source and field point, returns the lab time
/// \param fieldPoint The field point at the retarded time
/// \param ePosition The source position at the retarded time
/// \param tRet The retarded time being evaluated (in seconds)
/// \Returns The corresponding lab time
double CalcTimeFromRetardedTime(ROOT::Math::XYZPoint fieldPoint,
                                ROOT::Math::XYZPoint ePosition, double tRet);

/// Same function but using TVector3 apparatus
double CalcTimeFromRetardedTime(TVector3 fieldPoint, TVector3 ePosition,
                                double tRet);

/// Function for retrieving particle speed from kinetic energy
/// Basically just saves rewriting a few lines of code multiple times
/// \param T Particle kinetic energy in electronvolts
/// \param particleMass The particle mass in kilograms
/// \Returns The particle speed in metres per second
double GetSpeedFromKE(double T, double particleMass);

/// Function for calculating the gyroradius/Larmor radius/cyclotron radius of a
/// particle \param velocity Velocity vector of the particle in questions. Units
/// of metres per second. \param bField The magnetioc field vector at the
/// particle position in Tesla \param particleMass The particle mass in kg
double GetGyroradius(TVector3 velocity, TVector3 bField, double particleMass);

// Produces power spectrum with the desired normalisation
TGraph *MakePowerSpectrumNorm(const TGraph *grWave);

// Produdces power spectrum normalised to the mean squared amplitude of the time
// domain
TGraph *MakePowerSpectrumPeriodogram(const TGraph *grWave);

/// @brief Produces power spectrum
/// @param grWave Input time series graph
/// @return TGraph of periodogram
TGraph MakePowerSpectrumPeriodogram(const TGraph &grWave);

// Integrate the power spectrum
double IntegratePowerNorm(const TGraph *grFFT, Int_t firstBin = -1,
                          Int_t lastBin = -1);

// Implements a simple band pass filter
TGraph *BandPassFilter(const TGraph *grWave, const double minFreq,
                       const double maxFreq);

/// Function implementing a band pass filter on vectors of values
/// \param xVals A vector of equally spaced time values
/// \param yVals A vector of the corresponding y values
/// \param minFreq The value (in Hertz) below which to filter out frequencies
/// \param maxFreq The value (in Hertz) above which to filter out frequencies
/// \return A vector of the filtered y values
std::vector<double> BandPassFilter(std::vector<double> xVals,
                                   std::vector<double> yVals, double minFreq,
                                   double maxFreq);

/// @brief Given a set of 4 (x, y) values allow for interpolating at a value.
/// Use a 3rd order Lagrange interpolating polynomial
/// @param xVals A vector of 4 x values. These do not have to be evenly spaced
/// @param yVals A vector of the corresponding y values
/// @param xInterp The x value at which to interpolate
/// @return The interpolated y value
template <typename T>
T CubicInterpolation(std::vector<T> &xVals, std::vector<T> &yVals, T xInterp) {
  // Check we have the right number of values
  if (xVals.size() != 4 || yVals.size() != 4) {
    std::cout << "Invalid interpolation input size! Require 4 values.\n";
    return 0;
  }
  // Check that we have monotonically increasing x values
  for (unsigned int i{1}; i < xVals.size(); i++) {
    if (xVals.at(i) - xVals.at(i - 1) <= 0) {
      std::cout
          << "x values (for interpolation) do not increase monotonically!\n";
      return 0;
    }
  }

  // First calculate the Lagrange xInterpolating basis functions
  T l0{(xInterp - xVals[1]) * (xInterp - xVals[2]) * (xInterp - xVals[3]) /
       ((xVals[0] - xVals[1]) * (xVals[0] - xVals[2]) * (xVals[0] - xVals[3]))};
  T l1{(xInterp - xVals[0]) * (xInterp - xVals[2]) * (xInterp - xVals[3]) /
       ((xVals[1] - xVals[0]) * (xVals[1] - xVals[2]) * (xVals[1] - xVals[3]))};
  T l2{(xInterp - xVals[0]) * (xInterp - xVals[1]) * (xInterp - xVals[3]) /
       ((xVals[2] - xVals[0]) * (xVals[2] - xVals[1]) * (xVals[2] - xVals[3]))};
  T l3{(xInterp - xVals[0]) * (xInterp - xVals[1]) * (xInterp - xVals[2]) /
       ((xVals[3] - xVals[0]) * (xVals[3] - xVals[1]) * (xVals[3] - xVals[2]))};

  return l0 * yVals[0] + l1 * yVals[1] + l2 * yVals[2] + l3 * yVals[3];
}

/// @brief PDF of the Rayleigh distribution
/// @param x
/// @param sigma Scale parameter of the distribution
/// @return Probability density at x
double RayleighPDF(double x, double sigma);

/// @brief PDF of the Rayleigh distribution using long double
/// @param x
/// @param sigma Scale parameter of the distribution
/// @return Probability density at x
long double RayleighPDF(long double x, long double sigma);

/// @brief CDF for the Rayleigh distribution
/// @param x
/// @param sigma Scale parameter of the distribution
/// @return Cumulative distribution at x
double RayleighCDF(double x, double sigma);

/// @brief CDF for the Rayleigh distribution using long double
/// @param x
/// @param sigma Scale parameter of the distribution
/// @return Cumulative distribution at x
long double RayleighCDF(long double x, long double sigma);

double RayleighPDFFunc(double *x, double *par);

double RayleighCDFFunc(double *x, double *par);

void AddWhiteNoiseFrequencyDomainPowerNorm(TGraph *grIn, const double Teff,
                                           const int seed = 0);

/// \param BField A 3-vector describing the magnetic field at a point in space
/// \param charge Particle charge. Default is electron charge
/// \param energy Particle kinetic energy in electronvolts
/// \param mass Particle mass in kilograms. Default is electron rest mass
/// \Returns the cyclotron frequency
TVector3 calculate_omega(const TVector3 BField,
                         const double charge = -TMath::Qe(),
                         const double energy = 0.0, const double mass = ME);

/// Downmixes a time series of data with the appropriate frequency
/// \param grInput The input time series data
/// \param freq The frequency in Hertz with which to downmix
/// \Returns The downmixed time series
TGraph *DownmixInPhase(TGraph *grInput, const double freq);

/// Downmixes a time series of data with the appropriate frequency
/// \param grInput The input time series data
/// \param freq The frequency in Hertz with which to downmix
/// \Returns The downmixed time series
TGraph *DownmixQuadrature(TGraph *grInput, const double freq);

/// Scales all points of a TGraph by a constant
/// \param grInput The input graph
/// \param scale factor
/// \Returns grInput * scale
void ScaleGraph(TGraph *grInput, const double scale);

/// Calculates relativistic electron cyclotron frequency
/// \param KE The electron kinetic energy in eV
/// \param B Magnetic field strength in Tesla
double CalcCyclotronFreq(const double KE, const double B = 1.0);

/// Sums the points in a number of TGraphs
/// \param grInput A vector of the graphs to be summed
/// \Returns A graph containing the summed points
TGraph *SumGraphs(std::vector<TGraph *> grInput);

/// Downsamples a time series graph at a given sample rate using linear
/// interpolation \param grInput The input graph to be sampled \param sRate The
/// sample rate \Returns The downsampled graph
TGraph *SampleWaveform(TGraph *grInput, const double sRate);

/// Converts an input TGraph to a histogram
/// \param grInput Input graph to be converted
/// \return The converted histogram
TH1D *GraphToHistogram(TGraph *grInput);

/// Dowmixes, filters and downsamples a time series graph
/// \param grInput The inputted time domain voltage graph
/// \param downmixFreq The frequency at which to downmix, in Hertz
/// \param sampleRate The frequency at which to sample, in Hertz
/// \return The downmixed, filtered and sampled time domain graph
TGraph *SignalProcessGraph(TGraph *grInput, const double downmixFreq,
                           const double sampleRate);

/// Produces a graph of the FFT magnitudes
/// \param grInput The input real-valued time series data
/// \return The FFT magnitudes in frequency space
TGraph *MakeFFTMagGraph(TGraph *grInput);

/// Heaviside step function
/// \param x Input parameter
/// \return 0, for x <= 0 and 1 for x > 0
double HeavisideFunc(double x);

/// Rotates a given vector to the global coordinate system
/// \param v The vector to be rotated
/// \param xAx The X axis of the frame you are rotating from
/// \param yAx The Y axis of the frame you are rotating from
/// \param zAx The Z axis of the frame you are rotating from
/// \return The rotated vector
TVector3 RotateToGlobalCoords(TVector3 v, TVector3 xAx, TVector3 yAx,
                              TVector3 zAx);

/// Distribution for a skewed gaussian
/// \param x Point at which to evaluate the gaussian
/// \param A Scale factor
/// \param mu Centre of gaussian
/// \param sigma Width of gaussian
/// \param alpha Skewness
/// \return f(x) evaluated at x
double SkewedGaussian(double x, double A, double mu, double sigma,
                      double alpha);

/// @brief Function for a chirp signal
/// @param A Signal amplitude
/// @param t Time at which to generate signal [s]
/// @param phi0 Initial phase [radians]
/// @param f0 Initial frequency [Hz]
/// @param c Chirp rate [Hz s^-1]
/// @return The chirp function at the supplied time
double ChirpSignal(double A, double t, double phi0, double f0, double c);

/// Heaviside step function
/// \param x Input parameter
/// \Returns 1, for x > 0
/// \Returns 0, for x <= 0
double HeavisideFunc(double x);

/// Function for getting the zeros of the derivative of Bessel functions
/// \param n The order of the Bessel function being differentiated
/// \param m The zero of the derived function (must be > 0)
/// \Returns The specified root of the derivative of the nth Bessel function
double GetBesselPrimeZero(unsigned int n, unsigned int m);

/// @brief Larmor radiated power
/// @param ke Kinetic energy in eV
/// @param B Magnetic field strength in tesla
/// @param theta Pitch angle in radians
/// @param m Particle mass in kg
/// @return Larmor power in Watts
double CalcLarmorPower(double ke, double B, double theta, double m = ME);

/// @brief Rotate vector to different coordinate system
/// @param v Vector to be rotated
/// @param newX X axis of frame to be rotated to
/// @param newY Y axis of frame to be rotated to
/// @param newZ Z axis of frame to be rotated to
/// @return 3-vector of rotated vector
TVector3 RotateToCoords(TVector3 v, TVector3 newX, TVector3 newY,
                        TVector3 newZ);

/// @brief The linear sum of the power in a TGraph
/// @param gr A pointer to the input TGraph
/// @param firstBin The first bin to include in the sum
/// @param lastBin The last bin to include in the sum
/// @return The integral of the power
double SumPower(const TGraph *gr, int firstBin = -1, int lastBin = -1);

/// @brief The sum of the voltage squared in a waveform.
/// @param gr A pointer to the input TGraph
/// @param firstBin The first bin to include in the sum.
/// @param lastBin The last bin to include in the sum.
/// @return The value of the sum.
double SumVoltageSquared(const TGraph *gr, int firstBin = -1, int lastBin = -1);

/// @brief Computes the correlation of two arrays
/// @param length The length of the two arrays
/// @param oldY1 The first array in the correlation
/// @param oldY2 The second array in the correlation
/// @return The correlation as an array of *length* real numbers
double *GetCorrelation(int length, double *oldY1, double *oldY2);

/// @brief Computes the correlation of two TGraphs
/// @param gr1 The first graph in the correlation
/// @param gr2 The second graph in the correlation
/// @return A pointer to a tgraph containing the correaltion of gr1 and gr2
TGraph *GetCorrelationGraph(const TGraph *gr1, const TGraph *gr2,
                            int *zeroOffset = 0);

/// @brief Returns the normalised correlation of two TGraphs
/// @param gr1 The first graph in the correlation
/// @param gr2 The second graph in the correlation
/// @return A pointer to a TGraph containing the correlation of gr1 and
/// gr2 where each point is normalised by the number of valid samples in
/// the correlation and by the product of the RMS of the input graphs.
TGraph *GetNormalisedCorrelationGraph(const TGraph *gr1, const TGraph *gr2,
                                      int *zeroOffset = 0);

/// @brief Returns the normalised correlation of two TGraphs
/// @param gr1 The first input TGraph (must be zero meaned)
/// @param gr2 The second input TGraph (must be zero meaned)
/// @param zeroOffset A pointer to an integer where the sample corresponding to
/// zero offset will be stored
/// @param useDtRange A flag to enable the setting of a limited range of
/// deltat's to try
/// @param dtMin The minimum delta-t to include in the correlation, the maximum
/// delta-t to include in the correlation
/// @param dtMax
/// @return A pointer to a TGraph containing the correlation of <i>gr1</i> and
/// <i>gr2</i> where each point is normalised by the number of valid samples in
/// the correlation and by the product of the RMS of the input graphs.
TGraph *GetNormalisedCorrelationGraphTimeDomain(
    const TGraph *gr1, const TGraph *gr2, int *zeroOffset = 0,
    int useDtRange = 0, double dtMin = -1000, double dtMax = 1000);
}  // namespace rad

#endif
