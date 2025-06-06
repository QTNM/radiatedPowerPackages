// BasicFunctions.cxx

#include "BasicFunctions/BasicFunctions.h"

#include <cmath>
#include <iostream>
#include <vector>

#include "BasicFunctions/Constants.h"
#include "TAxis.h"
#include "TF1.h"
#include "TGraph.h"
#include "TH1.h"
#include "TH2.h"
#include "TMath.h"
#include "TRandom3.h"
#include "TVector3.h"

void rad::setGraphAttr(TGraph *gr) {
  gr->GetXaxis()->SetTitleSize(0.05);
  gr->GetYaxis()->SetTitleSize(0.05);
  gr->GetXaxis()->SetLabelSize(0.05);
  gr->GetYaxis()->SetLabelSize(0.05);
  gr->SetLineWidth(2);
}

void rad::setGraphAttr(std::unique_ptr<TGraph> &gr) {
  gr->GetXaxis()->SetTitleSize(0.05);
  gr->GetYaxis()->SetTitleSize(0.05);
  gr->GetXaxis()->SetLabelSize(0.05);
  gr->GetYaxis()->SetLabelSize(0.05);
  gr->SetLineWidth(2);
}

void rad::setGraphAttr(TGraph &gr) {
  gr.GetXaxis()->SetTitleSize(0.05);
  gr.GetYaxis()->SetTitleSize(0.05);
  gr.GetXaxis()->SetLabelSize(0.05);
  gr.GetYaxis()->SetLabelSize(0.05);
  gr.SetLineWidth(2);
  gr.SetLineColor(kBlack);
  gr.SetMarkerStyle(20);
  gr.SetMarkerColor(kBlack);
}

void rad::SetHistAttr(TH1 *h) {
  h->GetXaxis()->SetTitleSize(0.05);
  h->GetYaxis()->SetTitleSize(0.05);
  h->GetXaxis()->SetLabelSize(0.05);
  h->GetYaxis()->SetLabelSize(0.05);
  h->SetLineWidth(2);
}

void rad::SetHistAttr(TH1 &h) {
  h.GetXaxis()->SetTitleSize(0.05);
  h.GetYaxis()->SetTitleSize(0.05);
  h.GetXaxis()->SetLabelSize(0.05);
  h.GetYaxis()->SetLabelSize(0.05);
  h.SetLineWidth(2);
}

void rad::SetHistAttr(std::unique_ptr<TH1D> &h) {
  h->GetXaxis()->SetTitleSize(0.05);
  h->GetYaxis()->SetTitleSize(0.05);
  h->GetXaxis()->SetLabelSize(0.05);
  h->GetYaxis()->SetLabelSize(0.05);
  h->SetLineWidth(2);
}

void rad::SetHistAttr(TH2 *h) {
  h->GetXaxis()->SetTitleSize(0.05);
  h->GetYaxis()->SetTitleSize(0.05);
  h->GetZaxis()->SetTitleSize(0.05);
  h->GetXaxis()->SetLabelSize(0.05);
  h->GetYaxis()->SetLabelSize(0.05);
  h->GetZaxis()->SetLabelSize(0.05);
}

void rad::SetHistAttr(std::unique_ptr<TH2> &h) {
  h->GetXaxis()->SetTitleSize(0.05);
  h->GetYaxis()->SetTitleSize(0.05);
  h->GetZaxis()->SetTitleSize(0.05);
  h->GetXaxis()->SetLabelSize(0.05);
  h->GetYaxis()->SetLabelSize(0.05);
  h->GetZaxis()->SetLabelSize(0.05);
}

double rad::CalcAeHertzianDipole(const double wavelength,
                                 const ROOT::Math::XYZVector dipoleDir,
                                 const ROOT::Math::XYZPoint ePosition,
                                 const ROOT::Math::XYZPoint antennaPoint) {
  double Ae = 3 * pow(wavelength, 2) / (8 * TMath::Pi());
  const double psi =
      TMath::ACos(((ePosition - antennaPoint).Unit()).Dot(dipoleDir));
  Ae *= pow(TMath::Sin(psi), 2);
  return Ae;
}

double rad::CalcAlHertzianDipole(const double wavelength,
                                 const ROOT::Math::XYZVector dipoleDir,
                                 const ROOT::Math::XYZPoint ePosition,
                                 const ROOT::Math::XYZPoint antennaPoint) {
  double Al = wavelength * TMath::Sqrt(3 / (8 * TMath::Pi()));
  const double psi =
      TMath::ACos(((ePosition - antennaPoint).Unit()).Dot(dipoleDir));
  Al *= TMath::Sin(psi);
  return Al;
}

double rad::CalcRetardedTime(const ROOT::Math::XYZPoint fieldPoint,
                             const ROOT::Math::XYZPoint ePosition,
                             const double labTime) {
  double time =
      labTime - TMath::Sqrt((ePosition - fieldPoint).Mag2()) / TMath::C();
  return time;
}

double rad::CalcTimeFromRetardedTime(ROOT::Math::XYZPoint fieldPoint,
                                     ROOT::Math::XYZPoint ePosition,
                                     double tRet) {
  double time =
      tRet + TMath::Sqrt((ePosition - fieldPoint).Mag2()) / TMath::C();
  return time;
}

double rad::CalcTimeFromRetardedTime(TVector3 fieldPoint, TVector3 ePosition,
                                     double tRet) {
  double time = tRet + ((ePosition - fieldPoint).Mag() / TMath::C());
  return time;
}

double rad::GetSpeedFromKE(double T, double particleMass) {
  double gamma = T * TMath::Qe() / (ME * TMath::C() * TMath::C()) + 1;
  double betaSq = 1 - 1 / pow(gamma, 2);
  double speed = sqrt(betaSq) * TMath::C();
  return speed;
}

double rad::GetGyroradius(TVector3 velocity, TVector3 bField,
                          double particleMass) {
  double gamma{1 /
               sqrt(1 - velocity.Dot(velocity) / (TMath::C() * TMath::C()))};
  TVector3 vPerp{velocity - (velocity.Dot(bField.Unit()) * bField)};
  double rg{gamma * particleMass * vPerp.Mag() / (TMath::Qe() * bField.Mag())};
  return rg;
}

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

std::vector<double> rad::BandPassFilter(std::vector<double> xVals,
                                        std::vector<double> yVals,
                                        double minFreq, double maxFreq) {
  // Do some input checking
  if (xVals.size() != yVals.size()) {
    std::cout << "Your x and y values are not the same length. Are you sure "
                 "this is right?\n";
  }
  if (xVals.size() < 2 || yVals.size() < 2) {
    std::cout << "Invalid input value size. Returning standard output.\n";
    return std::vector<double>{0};
  } else {
    int length = yVals.size();
    double deltaT{xVals.at(1) - xVals.at(0)};
    FFTWComplex *theFFT = doFFT(length, &yVals[0]);

    int newLength{(length / 2) + 1};
    double deltaF{1.0 / (deltaT * length)};

    double tempF{0};
    for (int i{0}; i < newLength; i++) {
      if (tempF < minFreq || tempF > maxFreq) {
        theFFT[i].re = 0;
        theFFT[i].im = 0;
      }
      tempF += deltaT;
    }

    double *filteredVals = doInverseFFT(length, theFFT);
    std::vector<double> filteredValsVec(filteredVals, filteredVals + length);

    delete[] theFFT;
    delete[] filteredVals;
    return filteredValsVec;
  }
}

double rad::RayleighPDF(double x, double sigma) {
  double sigmaSq = sigma * sigma;
  double prob = (x / sigmaSq) * exp(-x * x / (2 * sigmaSq));
  return prob;
}

long double rad::RayleighPDF(long double x, long double sigma) {
  long double sigmaSq = sigma * sigma;
  long double prob = (x / sigmaSq) * exp(-x * x / (2 * sigmaSq));
  return prob;
}

double rad::RayleighCDF(double x, double sigma) {
  double f{1.0 - exp(-x * x / (2 * sigma * sigma))};
  return f;
}

long double rad::RayleighCDF(long double x, long double sigma) {
  long double f{1.0 - exp(-x * x / (2 * sigma * sigma))};
  return f;
}

double rad::RayleighPDFFunc(double *x, double *par) {
  double xx = x[0];
  double retVal = RayleighPDF(xx, par[0]);
  return retVal;
}

double rad::RayleighCDFFunc(double *x, double *par) {
  double xx = x[0];
  double retVal = RayleighCDF(xx, par[0]);
  return retVal;
}

void rad::AddWhiteNoiseFrequencyDomainPowerNorm(TGraph *grIn, const double Teff,
                                                const int seed) {
  gRandom->SetSeed(seed);
  const double sampleRate = 2 * grIn->GetPointX(grIn->GetN() - 1);
  const double deltaT = 1.0 / sampleRate;
  const double deltaF = grIn->GetPointX(1) - grIn->GetPointX(0);
  const double sigma = TMath::Sqrt(TMath::K() * Teff * sampleRate);
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

TVector3 rad::calculate_omega(const TVector3 BField, const double charge,
                              const double energy, const double mass) {
  double gamma_m0 = mass + energy * TMath::Qe() / pow(TMath::C(), 2);
  return (charge * BField * (1.0 / gamma_m0));
}

TGraph *rad::DownmixInPhase(TGraph *grInput, const double freq) {
  TGraph *grOut = new TGraph();
  setGraphAttr(grOut);
  for (int i = 0; i < grInput->GetN(); i++) {
    double time = grInput->GetPointX(i);
    grOut->SetPoint(
        i, time,
        grInput->GetPointY(i) * TMath::Cos(2 * TMath::Pi() * freq * time));
  }
  return grOut;
}

TGraph *rad::DownmixQuadrature(TGraph *grInput, const double freq) {
  TGraph *grOut = new TGraph();
  setGraphAttr(grOut);
  for (int i = 0; i < grInput->GetN(); i++) {
    double time = grInput->GetPointX(i);
    grOut->SetPoint(
        i, time,
        grInput->GetPointY(i) * TMath::Sin(2 * TMath::Pi() * freq * time));
  }
  return grOut;
}

void rad::ScaleGraph(TGraph *grInput, const double scale) {
  for (int i = 0; i < grInput->GetN(); i++) {
    grInput->SetPointY(i, grInput->GetPointY(i) * scale);
  }
}

double rad::CalcCyclotronFreq(const double KE, const double B) {
  double freq =
      TMath::Qe() * B / (ME + (KE * TMath::Qe() / pow(TMath::C(), 2)));
  return freq / (2 * TMath::Pi());
}

TGraph *rad::SumGraphs(std::vector<TGraph *> grInput) {
  TGraph *grOut = new TGraph();
  setGraphAttr(grOut);

  // First of all check that our graphs have the same x (time) spacing
  double testSpacing = grInput[0]->GetPointX(1) - grInput[0]->GetPointX(0);
  for (int iGr = 0; iGr < grInput.size(); iGr++) {
    double thisSpacing =
        grInput[iGr]->GetPointX(1) - grInput[iGr]->GetPointX(0);
    // Return empty graph if not
    if ((thisSpacing - testSpacing) / testSpacing > 1e-10) {
      std::cout << "Graphs do not have equivalent time spacing! -- returning "
                   "empty graph."
                << std::endl;
      return grOut;
    }
  }

  // The time series may be of different lengths
  // Need to determine the overlapping times between all the graphs
  double latestStart = -DBL_MAX;
  double earliestEnd = DBL_MAX;
  for (int iGr = 0; iGr < grInput.size(); iGr++) {
    if (grInput[iGr]->GetPointX(0) > latestStart)
      latestStart = grInput[iGr]->GetPointX(0);

    if (grInput[iGr]->GetPointX(grInput[iGr]->GetN() - 1) < earliestEnd)
      earliestEnd = grInput[iGr]->GetPointX(grInput[iGr]->GetN() - 1);
  }

  std::vector<int> startIndices;
  std::vector<int> endIndices;
  // Find the indices where these start and end times are satisfied
  for (int iGr = 0; iGr < grInput.size(); iGr++) {
    // Loop through points from start
    for (int iPnt = 0; iPnt < grInput[iGr]->GetN(); iPnt++) {
      if (grInput[iGr]->GetPointX(iPnt) == latestStart) {
        startIndices.push_back(iPnt);
        break;
      }
    }
    // Now loop through the points from the end to get the end
    for (int iPnt = grInput[iGr]->GetN() - 1; iPnt >= 0; iPnt--) {
      if (grInput[iGr]->GetPointX(iPnt) == earliestEnd) {
        endIndices.push_back(iPnt);
        break;
      }
    }  // Loop through points
  }

  // Now sum the graphs between the appropriate ranges
  int nPointsToSum = endIndices[0] - startIndices[0];
  for (int iPnt = 0; iPnt < nPointsToSum; iPnt++) {
    double xVal = grInput[0]->GetPointX(startIndices[0] + iPnt);
    double yVal = 0;
    for (int iGr = 0; iGr < grInput.size(); iGr++) {
      yVal += grInput[iGr]->GetPointY(startIndices[iGr] + iPnt);
    }
    grOut->SetPoint(grOut->GetN(), xVal, yVal);
  }

  return grOut;
}

TGraph *rad::SampleWaveform(TGraph *grInput, const double sRate) {
  TGraph *grOut = new TGraph();
  const double sampleSpacing = 1.0 / sRate;
  double sampleTime = grInput->GetPointX(0);

  for (int i = 0; i < grInput->GetN(); i++) {
    double time = grInput->GetPointX(i);
    if (time < sampleTime)
      continue;
    else if (i == 0) {
      double calcV = grInput->GetPointY(0);
      grOut->SetPoint(grOut->GetN(), sampleTime, calcV);
      sampleTime += sampleSpacing;
    } else {
      // Sample the distribution using linear interpolation
      double calcV = grInput->GetPointY(i - 1) +
                     (sampleTime - grInput->GetPointX(i - 1)) *
                         (grInput->GetPointY(i) - grInput->GetPointY(i - 1)) /
                         (time - grInput->GetPointX(i - 1));
      grOut->SetPoint(grOut->GetN(), sampleTime, calcV);
      sampleTime += sampleSpacing;
    }
  }

  return grOut;
}

TH1D *rad::GraphToHistogram(TGraph *grInput) {
  double binWidth = grInput->GetPointX(1) - grInput->GetPointX(0);
  double firstPoint = grInput->GetPointX(0);
  double lastPoint = grInput->GetPointX(grInput->GetN() - 1);
  TH1D *h = new TH1D("h", "", grInput->GetN(), firstPoint - binWidth / 2,
                     lastPoint + binWidth / 2);
  SetHistAttr(h);
  for (int i = 0; i < grInput->GetN(); i++) {
    h->SetBinContent(i + 1, grInput->GetPointY(i));
  }

  return h;
}

TGraph *rad::SignalProcessGraph(TGraph *grInput, const double downmixFreq,
                                const double sampleRate) {
  TGraph *grDM = DownmixInPhase(grInput, downmixFreq);
  TGraph *grS1 = SampleWaveform(grDM, 10 * sampleRate);
  delete grDM;
  TGraph *grF = BandPassFilter(grS1, 0.0, sampleRate / 2.0);
  delete grS1;
  TGraph *grS2 = SampleWaveform(grF, sampleRate);
  delete grF;

  return grS2;
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

double rad::HeavisideFunc(double x) {
  if (x > 0.0)
    return 1.0;
  else
    return 0.0;
}

double rad::GetBesselPrimeZero(unsigned int n, unsigned int m) {
  double zerosJ0Prime[5] = {3.8317, 7.0156, 10.1735, 13.3237, 16.4706};
  double zerosJ1Prime[5] = {1.8412, 5.3314, 8.5363, 11.7060, 14.8636};
  double zerosJ2Prime[5] = {3.0542, 6.7061, 9.9695, 13.1704, 16.3475};
  double zerosJ3Prime[5] = {4.2012, 8.0152, 11.3459, 14.5858, 17.7887};
  double zerosJ4Prime[5] = {5.3175, 9.2824, 12.6819, 15.9641, 19.1960};
  double zerosJ5Prime[5] = {6.4156, 10.5199, 13.9872, 17.3128, 20.5755};

  double p_prime_nm{0.0};
  if (m == 0) {
    std::cout
        << "Cannot have a zeroth zero of the function. Please choose m > 0."
        << std::endl;
    return p_prime_nm;
  } else if (n < 6) {
    if (n == 0) {
      p_prime_nm = zerosJ0Prime[m - 1];
    } else if (n == 1) {
      p_prime_nm = zerosJ1Prime[m - 1];
    } else if (n == 2) {
      p_prime_nm = zerosJ2Prime[m - 1];
    } else if (n == 3) {
      p_prime_nm = zerosJ3Prime[m - 1];
    } else if (n == 4) {
      p_prime_nm = zerosJ4Prime[m - 1];
    } else if (n == 5) {
      p_prime_nm = zerosJ5Prime[m - 1];
    }
    return p_prime_nm;
  } else {
    std::cout << "Currently don't have roots for this high n. Sorry!"
              << std::endl;
    return p_prime_nm;
  }
}

TVector3 rad::RotateToGlobalCoords(TVector3 v, TVector3 xAx, TVector3 yAx,
                                   TVector3 zAx) {
  // We are just transforming to the ROOT frame so our new axis coordinates
  // are just the unit vectors
  TVector3 unitX(1, 0, 0);
  TVector3 unitY(0, 1, 0);
  TVector3 unitZ(0, 0, 1);
  double x1Prime{unitX.Dot(xAx.Unit()) * v.X() + unitX.Dot(yAx.Unit()) * v.Y() +
                 unitX.Dot(zAx.Unit()) * v.Z()};
  double x2Prime{unitY.Dot(xAx.Unit()) * v.X() + unitY.Dot(yAx.Unit()) * v.Y() +
                 unitY.Dot(zAx.Unit()) * v.Z()};
  double x3Prime{unitZ.Dot(xAx.Unit()) * v.X() + unitZ.Dot(yAx.Unit()) * v.Y() +
                 unitZ.Dot(zAx.Unit()) * v.Z()};
  TVector3 newVector{x1Prime * unitX + x2Prime * unitY + x3Prime * unitZ};
  return newVector;
}

double rad::SkewedGaussian(double x, double A, double mu, double sigma,
                           double alpha) {
  // Evaluate gaussian
  double gaus{A * TMath::Gaus(x, mu, sigma)};
  // Evaluate
  double skew{1 + TMath::Erf(alpha * (x - mu) / (sigma * sqrt(2)))};
  return gaus * skew;
}

double rad::ChirpSignal(double A, double t, double phi0, double f0, double c) {
  return A * sin(phi0 + 2 * TMath::Pi() * (c * t * t / 2 + f0 * t));
}

TVector3 rad::RotateToCoords(TVector3 v, TVector3 newX, TVector3 newY,
                             TVector3 newZ) {
  // We are just transforming from the ROOT frame so our old axis coordinates
  // are just the unit vectors
  TVector3 oldX(1, 0, 0);
  TVector3 oldY(0, 1, 0);
  TVector3 oldZ(0, 0, 1);
  double pXPrime{newX.X() * v.X() + newY.X() * v.Y() + newZ.X() * v.Z()};
  double pYPrime{newX.Y() * v.X() + newY.Y() * v.Y() + newZ.Y() * v.Z()};
  double pZPrime{newX.Z() * v.X() + newY.Z() * v.Y() + newZ.Z() * v.Z()};
  TVector3 newVector(pXPrime, pYPrime, pZPrime);
  newVector = newVector.Unit();
  return newVector;
}

double rad::CalcLarmorPower(double ke, double B, double theta, double m) {
  const double f0{TMath::Qe() * B / (m * TMath::TwoPi())};
  const double beta{GetSpeedFromKE(ke, m) / TMath::C()};
  return TMath::TwoPi() * pow(TMath::Qe() * f0 * beta * sin(theta), 2) /
         ((3 * EPSILON0 * TMath::C()) * (1 - beta * beta));
}

double rad::SumPower(const TGraph *gr, int firstBin, int lastBin) {
  double integral{0};
  double deltaF = gr->GetPointX(1) - gr->GetPointX(0);
  if (firstBin < 0) firstBin = 0;
  if (lastBin < 0) lastBin = gr->GetN() - 1;
  for (int i = firstBin; i <= lastBin; i++) {
    integral += gr->GetPointY(i);
  }
  return integral;
}

double rad::SumVoltageSquared(const TGraph *gr, int firstBin, int lastBin) {
  double integral{0};
  if (firstBin < 0) firstBin = 0;
  if (lastBin < 0) lastBin = gr->GetN() - 1;
  for (int i = firstBin; i <= lastBin; i++) {
    integral += gr->GetPointY(i) * gr->GetPointY(i);
  }
  return integral;
}

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