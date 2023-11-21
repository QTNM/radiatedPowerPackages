/*
  ButterworthFilter.h

  Class describing Butterworth filters up to order 10

  S. Jones
  20/11/2023
*/

#include <boost/math/special_functions/binomial.hpp>
#include <boost/numeric/ublas/io.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/vector.hpp>
#include <cmath>
#include <iostream>
#include <vector>

#ifndef BUTTERWORTH_FILTER_H
#define BUTTERWORTH_FILTER_H

namespace rad {

typedef boost::numeric::ublas::vector<double, std::vector<double>> bvector;
typedef boost::numeric::ublas::matrix<double, boost::numeric::ublas::row_major,
                                      std::vector<double>>
    bmatrix;

class ButterworthFilter {
 public:
  /// @brief Parametrised constructor
  /// @param order Filter order
  /// @param cutoffFreq Cut-off frequency in Hertz
  /// @param sampleRate Sample rate in Hertz
  ButterworthFilter(unsigned int order, double cutoffFreq, double sampleRate);

  /// @brief Getter function for filter order
  /// @return Integer filter order
  unsigned int GetOrder() { return n; }

  /// @brief Getter function for cut-off frequency
  /// @return Desired cutoff frequency in Hertz
  double GetCutoffFrequency() { return fc; }

  /// @brief Getter function for transfer function coeffs
  /// @return Vector of transfer function numerator coeffs
  std::vector<double> GetaVec() { return a; }

  /// @brief Getter function for transfer function coeffs
  /// @return Vector of transfer function denominator coeffs
  std::vector<double> GetbVec() { return b; }

 private:
  unsigned int n;  // Filter order
  double fc;       // Cutoff frequency in Hertz
  double fs;       // Sample rate in Hertz

  std::vector<double> a;  // Vector of coefficients
  std::vector<double> b;  // Vector of coefficients

  std::vector<double> xk;  // Vector of coefficients
  std::vector<double> yk;  // Vector of coefficients

  /// @brief Gets coefficients from analog transfer function
  /// @return Vector of dimension n + 1
  bvector GetBCoeffs();

  /// @brief
  /// @return Matrix of size (n + 1; n + 1)
  bmatrix GetPMatrix();

  /// @brief Calculates matrix from pre-warping values
  /// @param c Parameter related to pre-warping of frequencies
  /// @return Diagonal matrix of size (n + 1; n + 1)
  bmatrix GetTcMatrix(double c);

  /// @brief Set the a vector
  void SetaVec();

  /// @brief Set the b vector
  void SetbVec();

  /// @brief Sets difference equation coeffs
  void SetDifferenceEqnCoeffs();
};
}  // namespace rad

#endif