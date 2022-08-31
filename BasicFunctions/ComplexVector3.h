/*
  ComplexVector3.h

  Class supporting the use of complex three vector
*/

#ifndef COMPLEX_VECTOR_THREE_H
#define COMPLEX_VECTOR_THREE_H

#include <complex>

namespace rad {

  class ComplexVector3 {

  private:
    // The components
    std::complex<double> mX, mY, mZ;

  public:
    /// Parametrised constructor taking three complex numbers
    ComplexVector3(std::complex<double> x, std::complex<double> y, std::complex<double> z) : mX{x}, mY{y}, mZ{z} {}
  };
}

#endif
