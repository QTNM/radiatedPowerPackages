/*
  ButterworthFilter.cxx
*/

#include "BasicFunctions/ButterworthFilter.h"

rad::ButterworthFilter::ButterworthFilter(unsigned int order, double cutoffFreq,
                                          double sampleRate)
    : n(order), fc(cutoffFreq), fs(sampleRate) {
  // Do the various calculations
  SetaVec();
  SetbVec();
  SetDifferenceEqnCoeffs();
}

void rad::ButterworthFilter::SetaVec() {
  a.resize(n + 1);
  for (size_t i{0}; i < a.size(); i++) {
    a.at(i) = boost::math::binomial_coefficient<double>(n, i);
  }
}

void rad::ButterworthFilter::SetbVec() {
  bvector B{GetBCoeffs()};
  const double theta_c{2 * M_PI * fc / fs};
  const double c{1.0 / tan(theta_c / 2)};
  bmatrix P{GetPMatrix()};
  bmatrix Tc{GetTcMatrix(c)};

  // Multiply matrices together to get the correct options
  bmatrix interMatrix{boost::numeric::ublas::prod(Tc, P)};
  bvector bVec{boost::numeric::ublas::prod(B, interMatrix)};

  b.resize(n + 1);
  for (size_t i{0}; i < b.size(); i++) {
    b.at(i) = bVec(i);
  }
}

rad::bvector rad::ButterworthFilter::GetBCoeffs() {
  std::vector<double> bVec;
  if (n == 1) {
    bVec = {1, 1};
  } else if (n == 2) {
    bVec = {1, sqrt(2), 1};
  } else if (n == 3) {
    bVec = {1, 2, 2, 1};
  } else if (n == 4) {
    bVec = {1, 2.6131, 3.4142, 2.6131, 1};
  } else if (n == 5) {
    bVec = {1, 3.2361, 5.2361, 5.2361, 3.2361, 1};
  } else if (n == 6) {
    bVec = {1, 3.8637, 7.4641, 9.1416, 7.4641, 3.8637, 1};
  } else if (n == 7) {
    bVec = {1, 4.4940, 10.0978, 14.5918, 14.5918, 10.0978, 4.4940, 1};
  } else if (n == 8) {
    bVec = {1, 5.1258, 13.1371, 21.8462, 25.6884, 21.8462, 5.1258, 13.1371, 1};
  } else if (n == 9) {
    bVec = {1,       5.7588,  16.5817, 31.1634, 41.9864,
            41.9864, 31.1634, 16.5817, 5.7588,  1};
  } else if (n == 10) {
    bVec = {1,       6.3925,  20.4317, 42.8021, 64.8824, 74.2334,
            64.8824, 42.8021, 20.4317, 6.3925,  1};
  }

  bvector B(bVec);
  return B;
}

rad::bmatrix rad::ButterworthFilter::GetPMatrix() {
  bmatrix P(n + 1, n + 1);
  // First set all elements in first column to zero
  for (unsigned int i{0}; i < P.size1(); ++i) {
    P(i, 0) = 1;
  }

  // Then set the first row of elements
  for (unsigned int j{0}; j < P.size2(); ++j) {
    P(0, j) = boost::math::binomial_coefficient<double>(n, j);
  }

  // Now do the remaining elements of the matrix
  for (unsigned int i{1}; i < P.size1(); ++i) {
    for (unsigned int j{1}; j < P.size2(); ++j) {
      P(i, j) = P(i - 1, j) - P(i, j - 1) - P(i - 1, j - 1);
    }
  }

  return P;
}

rad::bmatrix rad::ButterworthFilter::GetTcMatrix(double c) {
  bmatrix Tc(n + 1, n + 1);
  int power{0};
  for (unsigned int i{0}; i < Tc.size1(); ++i) {
    for (unsigned int j{0}; j < Tc.size2(); ++j) {
      if (i == j) {
        Tc(i, j) = pow(c, power);
        power++;
      } else
        Tc(i, j) = 0;
    }
  }
  return Tc;
}

void rad::ButterworthFilter::SetDifferenceEqnCoeffs() {
  xk.resize(n + 1, 0);
  yk.resize(n + 1, 0);
  for (int i{0}; i < xk.size(); i++) {
    if (i == 0) {
      xk.at(i) = a.at(i) / b.at(0);
    } else {
      xk.at(i) = a.at(i) / b.at(0);
      yk.at(i) = -b.at(i) / b.at(0);
    }
  }
}