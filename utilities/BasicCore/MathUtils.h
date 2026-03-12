// MathUtils.h
// Simple mathematical utility functions - header-only with no external
// dependencies

#ifndef MATHUTILS_H
#define MATHUTILS_H

#include <iostream>
#include <vector>

#include "Constants.h"

namespace rad {

/// @brief Convert degrees to radians
/// @param degrees Angle in degrees
/// @return Angle in radians
inline constexpr double DegToRad(double degrees) {
  return degrees * PI / 180.0;
}

/// @brief Convert radians to degrees
/// @param radians Angle in radians
/// @return Angle in degrees
inline constexpr double RadToDeg(double radians) {
  return radians * 180.0 / PI;
}

/// @brief Given 4 (x,y) values, interpolate using 3rd order Lagrange polynomial
/// @param xVals Vector of 4 x values (not necessarily evenly spaced)
/// @param yVals Vector of 4 corresponding y values
/// @param xInterp The x value at which to interpolate
/// @return Interpolated y value
template <typename T>
T CubicInterpolation(std::vector<T> &xVals, std::vector<T> &yVals, T xInterp) {
  // Check we have the right number of values
  if (xVals.size() != 4 || yVals.size() != 4) {
    std::cout << "Invalid interpolation input size! Require 4 values.\\n";
    return 0;
  }
  // Check that we have monotonically increasing x values
  for (unsigned int i{1}; i < xVals.size(); i++) {
    if (xVals.at(i) - xVals.at(i - 1) <= 0) {
      std::cout
          << "x values (for interpolation) do not increase monotonically!\\n";
      return 0;
    }
  }

  // Calculate the Lagrange interpolating basis functions
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
/// @param x location at which to calculate PDF
/// @param sigma Scale parameter of the distribution
/// @return Probability density at x
inline double RayleighPDF(double x, double sigma) {
  return (x / (sigma * sigma)) * exp(-x * x / (2 * (sigma * sigma)));
}

/// @brief PDF of the Rayleigh distribution
/// @param x Location at which to calculate PDF
/// @param sigma Scale parameter of the distribution
/// @return Probability density at x
inline long double RayleighPDF(long double x, long double sigma) {
  return (x / (sigma * sigma)) * exp(-x * x / (2 * (sigma * sigma)));
}

/// @brief CDF for the Rayleigh distribution
/// @param x Location at which to calculate PDF
/// @param sigma Scale parameter of the distribution
/// @return Cumulative distribution at x
inline double RayleighCDF(double x, double sigma) {
  return 1.0 - exp(-x * x / (2 * sigma * sigma));
}

/// @brief CDF for the Rayleigh distribution using long double
/// @param x Location at which to calculate PDF
/// @param sigma Scale parameter of the distribution
/// @return Cumulative distribution at x
inline long double RayleighCDF(long double x, long double sigma) {
  return 1.0 - exp(-x * x / (2 * sigma * sigma));
}

/// @brief Rayleigh PDF function wrapper for ROOT TF1
/// @param x Input variable array
/// @param par Parameter array
/// @return PDF value
inline double RayleighPDFFunc(double *x, double *par) {
  return RayleighPDF(x[0], par[0]);
}

/// @brief Rayleigh CDF function wrapper for ROOT TF1
/// @param x Input variable array
/// @param par Parameter array
/// @return CDF value
inline double RayleighCDFFunc(double *x, double *par) {
  return RayleighCDF(x[0], par[0]);
}

/// @brief Heaviside step function
/// @param x Input parameter
/// @return 0 for x < 0, 1 for x >= 0
inline double HeavisideFunc(double x) {
  if (x >= 0.0)
    return 1.0;
  else
    return 0.0;
}

/// @brief Get zeros of the derivative of Bessel functions
/// @param n Order of the Bessel function being differentiated
/// @param m The zero of the derived function (must be > 0)
/// @return The specified root of the derivative of the nth Bessel function
inline double GetBesselPrimeZero(unsigned int n, unsigned int m) {
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

/// @brief Distribution for a skewed gaussian
/// @param x Point at which to evaluate the gaussian
/// @param A Scale factor
/// @param mu Centre of gaussian
/// @param sigma Width of gaussian
/// @param alpha Skewness
/// @return f(x) evaluated at x
inline double SkewedGaussian(double x, double A, double mu, double sigma,
                             double alpha) {
  // Evaluate gaussian - use standard math instead of TMath
  double gaus{A * exp(-0.5 * pow((x - mu) / sigma, 2)) /
              (sigma * sqrt(2 * PI))};
  // Evaluate error function approximation
  double skew{1 + erf(alpha * (x - mu) / (sigma * sqrt(2)))};
  return gaus * skew;
}

}  // namespace rad

#endif