/// CorrelationAnalysis.cxx - ROOT-based correlation analysis functions
#include "ROOTUtils/CorrelationAnalysis.h"

#include "SignalUtils/FourierTransforms.h"
#include "SignalUtils/FFTWComplex.h" 
#include "TMath.h"

double *rad::GetCorrelation(int length, double *oldY1, double *oldY2) {
  double *newY1 = new double[length];
  FFTWComplex *theFFT1 = doFFT(length, oldY1);
  FFTWComplex *theFFT2 = doFFT(length, oldY2);
  int newLength = (length / 2) + 1;
  FFTWComplex *tempStep = new FFTWComplex[newLength];
  int no2 = length >> 1;
  for (int i{0}; i < newLength; i++) {
    double reFFT1{theFFT1[i].re};
    double imFFT1{theFFT1[i].im};
    double reFFT2{theFFT2[i].re};
    double imFFT2{theFFT2[i].im};
    tempStep[i].re = (reFFT1 * reFFT2 + imFFT1 * imFFT2) / double(no2 / 2);
    tempStep[i].re = (imFFT1 * reFFT2 - reFFT1 * imFFT2) / double(no2 / 2);
  }

  double *theOutput = doInverseFFT(length, tempStep);
  delete[] theFFT1;
  delete[] theFFT2;
  delete[] tempStep;
  return theOutput;
}

TGraph *rad::GetCorrelationGraph(const TGraph *gr1, const TGraph *gr2,
                                 int *zeroOffset) {
  // Now we'll extend this up to a power of 2
  int length = gr1->GetN();
  int length2 = gr2->GetN();

  int N = int(TMath::Power(2, int(TMath::Log2(length)) + 2));
  if (N < length2) N = int(TMath::Power(2, int(TMath::Log2(length2)) + 2));

  // Will really assume that N's are equal for now
  int firstRealSamp = (N - length) / 2;

  double *oldY1 = new double[N];
  double *oldY2 = new double[N];

  double x, y;
  double x2, y2;
  gr1->GetPoint(1, x2, y2);
  gr1->GetPoint(0, x, y);
  double deltaT = x2 - x;
  double firstX = x;

  gr2->GetPoint(0, x2, y2);
  double waveOffset = firstX - x2;

  for (int i = 0; i < N; i++) {
    if (i < firstRealSamp || i >= firstRealSamp + length)
      y = 0;
    else {
      gr1->GetPoint(i - firstRealSamp, x, y);
    }
    oldY1[i] = y;

    if (i < firstRealSamp || i >= firstRealSamp + length2)
      y = 0;
    else {
      gr2->GetPoint(i - firstRealSamp, x, y);
    }
    oldY2[i] = y;
  }

  if (zeroOffset) {
    *zeroOffset = N / 2;
    (*zeroOffset) += Int_t(waveOffset / deltaT);
  }

  double *xVals = new double[N];
  double *yVals = new double[N];
  double *corVals = GetCorrelation(N, oldY1, oldY2);
  for (int i = 0; i < N; i++) {
    if (i < N / 2) {
      // Positive
      xVals[i + (N / 2)] = (i * deltaT) + waveOffset;
      yVals[i + (N / 2)] = corVals[i];
    } else {
      // Negative
      xVals[i - (N / 2)] = ((i - N) * deltaT) + waveOffset;
      yVals[i - (N / 2)] = corVals[i];
    }
  }

  TGraph *grCor = new TGraph(N, xVals, yVals);
  delete[] oldY1;
  delete[] oldY2;
  delete[] xVals;
  delete[] yVals;
  delete[] corVals;

  return grCor;
}

TGraph *rad::GetNormalisedCorrelationGraph(const TGraph *gr1, const TGraph *gr2,
                                           int *zeroOffset) {
  // Will also assume these graphs are zero meaned
  // Now we'll extend this up to a power of 2
  int length{gr1->GetN()};
  Double_t *y1{gr1->GetY()};
  int length2{gr2->GetN()};
  Double_t *y2{gr2->GetY()};
  Double_t denom{gr1->GetRMS(2) * gr2->GetRMS(2)};

  int N = int(TMath::Power(2, int(TMath::Log2(length)) + 2));
  if (N < length2) N = int(TMath::Power(2, int(TMath::Log2(length2)) + 2));

  // Will really assume that N's are equal for now
  int firstRealSamp = 1 + (N - 2 * length) / 2;
  int lastRealSamp = firstRealSamp + 2 * (length - 1);
  TGraph *grCor = GetCorrelationGraph(gr1, gr2, zeroOffset);
  Double_t *corVal = grCor->GetY();
  Double_t norm1{0};
  Double_t norm2{0};

  for (int i{0}; i < N; i++) {
    if (i >= firstRealSamp && i <= lastRealSamp) {
      if (i <= N / 2) {
        norm1 += (y1[i - firstRealSamp] * y1[i - firstRealSamp]);
        norm2 += (y2[length - 1 - (i - firstRealSamp)] *
                  y2[length - 1 - (i - firstRealSamp)]);
        int effN = 1 + (i - firstRealSamp);
        corVal[i] /= (sqrt(effN) * denom);
      } else if (i < N - 1) {
        norm1 -= (y1[i - 1 - (N / 2)] * y1[i - 1 - (N / 2)]);
        norm2 -= (y2[length - (i - (N / 2))] * y2[length - (i - (N / 2))]);
        int effN = (1 + lastRealSamp - i);
        corVal[i] /= (sqrt(effN) * denom);
      }
    }
  }

  return grCor;
}

TGraph *rad::GetNormalisedCorrelationGraphTimeDomain(
    const TGraph *gr1, const TGraph *gr2, int *zeroOffset, int useDtRange,
    double dtMin, double dtMax) {
  // Will also assume these graphs are zero meaned... may fix this assumption
  // Now we'll extend this up to a power of 2
  int length{gr1->GetN()};
  double *y1{gr1->GetY()};
  int length2{gr2->GetN()};
  if (length2 < length) length = length2;
  double *y2{gr2->GetY()};
  double denom{gr1->GetRMS(2) * gr2->GetRMS(2)};

  double *x1{gr1->GetX()};
  double *x2{gr2->GetX()};

  double deltaT = x1[1] - x1[0];
  double waveOffset = x1[0] - x2[0];

  int N{2 * length - 1};

  if (zeroOffset) {
    *zeroOffset = N / 2;
    (*zeroOffset) += Int_t(waveOffset / deltaT);
  }

  // Will really assume that N's are equal for now
  int firstRealSamp = 0;
  int lastRealSamp = N - 1;
  int minDtIndex = 0;
  int maxDtIndex = N - 1;
  if (useDtRange) {
    minDtIndex = TMath::Floor((dtMin - waveOffset) / deltaT) + (N / 2);
    if (minDtIndex < 0) minDtIndex = 0;
    maxDtIndex = TMath::Ceil((dtMax - waveOffset) / deltaT) + (N / 2);
    if (maxDtIndex < 0) maxDtIndex = 0;
  }

  double *xVals = new double[N];
  double *corVals = new double[N];
  for (int i = minDtIndex; i <= maxDtIndex; i++) {
    int dtIndex = (i - minDtIndex);

    xVals[dtIndex] = ((i - N / 2) * deltaT) + waveOffset;
    corVals[dtIndex] = 0;
    if (i >= firstRealSamp && i <= lastRealSamp) {
      Int_t firstIndex = (i - firstRealSamp);
      Int_t secondIndex = length - 1;
      if (firstIndex > length - 1) {
        int offset = firstIndex - (length - 1);
        firstIndex = length - 1;
        secondIndex -= offset;
      }

      Int_t numSamples = 0;
      for (; firstIndex >= 0 && secondIndex >= 0; firstIndex--) {
        corVals[dtIndex] += y1[firstIndex] * y2[secondIndex];
        numSamples++;
        secondIndex--;
      }
      corVals[dtIndex] /= denom * sqrt(numSamples);
    }
  }

  TGraph *grCor = new TGraph((maxDtIndex - minDtIndex) + 1, xVals, corVals);
  delete[] xVals;
  delete[] corVals;
  return grCor;
}