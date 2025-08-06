/// PureFunctions.cxx - Pure mathematical functions without ROOT dependencies
#include "utilities/BasicFunctions/BasicFunctions.h"

#include <cmath>
#include <iostream>
#include <vector>

#include "utilities/BasicCore/Constants.h"
#include "utilities/SignalUtils/FourierTransforms.h"
#include "utilities/SignalUtils/FFTWComplex.h"

/// Helper function for GetSpeedFromKE
double GetSpeedFromKE(double T, double particleMass) {
  double gamma = T * rad::QE / (particleMass * rad::C * rad::C) + 1;
  double betaSq = 1 - 1 / (gamma * gamma);
  double speed = sqrt(betaSq) * rad::C;
  return speed;
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
      tempF += deltaF;  // Fixed: was deltaT, should be deltaF
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

double rad::SkewedGaussian(double x, double A, double mu, double sigma,
                           double alpha) {
  // Evaluate gaussian - use standard math instead of TMath
  double gaus = A * exp(-0.5 * pow((x - mu) / sigma, 2)) / (sigma * sqrt(2 * PI));
  // Evaluate error function approximation
  double skew = 1 + erf(alpha * (x - mu) / (sigma * sqrt(2)));
  return gaus * skew;
}

double rad::ChirpSignal(double A, double t, double phi0, double f0, double c) {
  return A * sin(phi0 + 2 * PI * (c * t * t / 2 + f0 * t));
}

double rad::CalcLarmorPower(double ke, double B, double theta, double m) {
  const double f0 = QE * B / (m * (2 * PI));
  const double beta = GetSpeedFromKE(ke, m) / C;
  return (2 * PI) * pow(QE * f0 * beta * sin(theta), 2) /
         ((3 * EPSILON0 * C) * (1 - beta * beta));
}