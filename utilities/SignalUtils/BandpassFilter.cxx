#include "utilities/SignalUtils/BandpassFilter.h"

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