/// FFTAnalysis.cxx - ROOT-based FFT and frequency domain analysis
#include "utilities/ROOTUtils/FFTAnalysis.h"

#include "utilities/BasicCore/Constants.h"
#include "utilities/ROOTUtils/GraphUtils.h"
#include "utilities/SignalUtils/FourierTransforms.h"
#include "utilities/SignalUtils/FFTWComplex.h"
#include <cmath>
#include "TF1.h"
#include "TRandom3.h"

using namespace rad;

// Very similar to the FFTtools implementation but without the scaling of the x
// axis the MHz
TGraph *rad::MakePowerSpectrumNorm(const TGraph *grWave) {
  double *oldY = grWave->GetY();
  double *oldX = grWave->GetX();
  double deltaT = oldX[1] - oldX[0];
  int length = grWave->GetN();
  FFTWComplex *theFFT = doFFT(length, oldY);

  int newLength = (length / 2) + 1;
  double *newY = new double[newLength];
  double *newX = new double[newLength];

  double deltaF = 1 / (deltaT * length);

  double tempF = 0;
  for (int i = 0; i < newLength; i++) {
    float power = pow(getAbs(theFFT[i]), 2);
    if (i > 0 && i < newLength - 1) power *= 2;  // account for symmetry
    power *= deltaT / (length);  // For time-integral squared amplitude
    power /= deltaF;             // Just to normalise bin-widths
    // Ends up the same as dt^2, need to integrate the power (multiply by df)
    // to get a meaningful number out.
    newX[i] = tempF;
    newY[i] = power;
    tempF += deltaF;
  }

  TGraph *grPower = new TGraph(newLength, newX, newY);
  setGraphAttr(grPower);
  delete[] theFFT;
  delete[] newY;
  delete[] newX;
  return grPower;
}

TGraph *rad::MakePowerSpectrumPeriodogram(const TGraph *grWave) {
  double *oldY = grWave->GetY();
  double *oldX = grWave->GetX();
  double deltaT = oldX[1] - oldX[0];
  int length = grWave->GetN();
  FFTWComplex *theFFT = doFFT(length, oldY);
  double lengthDub = (double)length;
  int newLength = (length / 2) + 1;
  double *newY = new double[newLength];
  double *newX = new double[newLength];

  double deltaF = 1 / (deltaT * length);

  double tempF = 0;
  for (int i = 0; i < newLength; i++) {
    float power = pow(getAbs(theFFT[i]), 2);
    if (i > 0 && i < newLength - 1) power *= 2;  // account for symmetry
    double scale = lengthDub * lengthDub;
    power /= scale;
    newX[i] = tempF;
    newY[i] = power;
    tempF += deltaF;
  }

  TGraph *grPower = new TGraph(newLength, newX, newY);
  setGraphAttr(grPower);
  grPower->GetXaxis()->SetTitle("Frequency [Hz]");
  delete[] theFFT;
  delete[] newY;
  delete[] newX;
  return grPower;
}

TGraph rad::MakePowerSpectrumPeriodogram(const TGraph &grWave) {
  double *oldY = grWave.GetY();
  double *oldX = grWave.GetX();
  double deltaT = oldX[1] - oldX[0];
  int length = grWave.GetN();
  FFTWComplex *theFFT = doFFT(length, oldY);
  double lengthDub = (double)length;
  int newLength = (length / 2) + 1;
  double *newY = new double[newLength];
  double *newX = new double[newLength];

  double deltaF = 1 / (deltaT * length);

  double tempF = 0;
  for (int i = 0; i < newLength; i++) {
    float power = pow(getAbs(theFFT[i]), 2);
    if (i > 0 && i < newLength - 1) power *= 2;  // account for symmetry
    double scale = lengthDub * lengthDub;
    power /= scale;
    newX[i] = tempF;
    newY[i] = power;
    tempF += deltaF;
  }

  TGraph grPower(newLength, newX, newY);
  setGraphAttr(grPower);
  grPower.GetXaxis()->SetTitle("Frequency [Hz]");
  delete[] theFFT;
  delete[] newY;
  delete[] newX;
  return grPower;
}

double rad::IntegratePowerNorm(const TGraph *grFFT, Int_t firstBin,
                               Int_t lastBin) {
  double integral{0};
  double freq{}, power{};
  // Multiply by frequency bin width
  double deltaF = grFFT->GetPointX(1) - grFFT->GetPointX(0);
  for (int i = firstBin; i <= lastBin; i++) {
    integral += grFFT->GetPointY(i) * deltaF;
  }
  integral *= deltaF;
  return integral;
}

// Re-implementation of the filter from FFTtools but without the conversion
// factors
TGraph *rad::BandPassFilter(const TGraph *grWave, const double minFreq,
                            const double maxFreq) {
  double *oldY = grWave->GetY();
  double *oldX = grWave->GetX();
  double deltaT = oldX[1] - oldX[0];
  int length = grWave->GetN();
  FFTWComplex *theFFT = doFFT(length, oldY);

  int newLength = (length / 2) + 1;
  double deltaF = 1 / (deltaT * length);  // Hz

  double tempF = 0;
  for (int i = 0; i < newLength; i++) {
    if (tempF < minFreq || tempF > maxFreq) {
      theFFT[i].re = 0;
      theFFT[i].im = 0;
    }
    tempF += deltaF;
  }

  double *filteredVals = doInverseFFT(length, theFFT);

  TGraph *grFiltered = new TGraph(length, oldX, filteredVals);
  delete[] theFFT;
  delete[] filteredVals;
  return grFiltered;
}

TGraph *rad::MakeFFTMagGraph(TGraph *grInput) {
  double *oldY = grInput->GetY();
  double *oldX = grInput->GetX();
  double deltaT = oldX[1] - oldX[0];
  int length = grInput->GetN();
  FFTWComplex *theFFT = doFFT(length, oldY);
  double lengthDub = (double)length;
  int newLength = (length / 2) + 1;
  double *newY = new double[newLength];
  double *newX = new double[newLength];

  double deltaF = 1 / (deltaT * length);

  double tempF = 0;
  for (int i = 0; i < newLength; i++) {
    float mag = getAbs(theFFT[i]);
    newX[i] = tempF;
    newY[i] = mag;
    tempF += deltaF;
  }

  TGraph *grMag = new TGraph(newLength, newX, newY);
  setGraphAttr(grMag);

  delete[] theFFT;
  delete[] newY;
  delete[] newX;
  return grMag;
}

/// PDF for Rayleigh distribution
double RayleighPDFFunc(double *x, double *par) {
  double xx = x[0];
  double sigma = par[0];
  double sigmaSq = sigma * sigma;
  double prob = (xx / sigmaSq) * exp(-xx * xx / (2 * sigmaSq));
  return prob;
}

void rad::AddWhiteNoiseFrequencyDomainPowerNorm(TGraph *grIn, const double Teff,
                                                const int seed) {
  gRandom->SetSeed(seed);
  const double sampleRate = 2 * grIn->GetPointX(grIn->GetN() - 1);
  const double deltaT = 1.0 / sampleRate;
  const double deltaF = grIn->GetPointX(1) - grIn->GetPointX(0);
  const double sigma = std::sqrt(K_B * Teff * sampleRate);
  TF1 *f1 = new TF1("f1", RayleighPDFFunc, 0, 4 * sigma, 1);
  f1->SetParameter(0, sigma);

  // Now loop through the bins of the graph and add the noise
  for (int i = 0; i < grIn->GetN(); i++) {
    double noise =
        pow(f1->GetRandom() * sqrt(0.5), 2) * (1 / deltaF) * (deltaT);
    grIn->SetPointY(i, grIn->GetPointY(i) + noise);
  }

  delete f1;
}