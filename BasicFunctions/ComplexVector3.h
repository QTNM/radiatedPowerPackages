/*
  ComplexVector3.h

  Class supporting the use of complex three vector
*/

#ifndef COMPLEX_VECTOR_THREE_H
#define COMPLEX_VECTOR_THREE_H

#include <complex>

#include "TVector3.h"

namespace rad {

  class ComplexVector3 {

  private:
    // The components
    std::complex<double> mX, mY, mZ;

  public:
    /// Parametrised constructor taking three complex numbers
    ComplexVector3(std::complex<double> x, std::complex<double> y, std::complex<double> z) : mX{x}, mY{y}, mZ{z} {}

    /// Parametrised constructor from a real three vector
    ComplexVector3(TVector3 vec);
    
    /// Functions to return individual components
    std::complex<double> X() { return mX; }
    std::complex<double> Y() { return mY; }
    std::complex<double> Z() { return mZ; }

    /// Functions to set components after creation
    void SetX(std::complex<double> x) { mX = x; }
    void SetY(std::complex<double> y) { mY = y; }
    void SetZ(std::complex<double> z) { mZ = z; }

    /// Scalar product
    /// \param vec The vector with which to compute the scalar product
    /// \Returns The complex scalar product of the two vectors
    std::complex<double> Dot(ComplexVector3 vec);

    /// Vector product
    /// \param vec The vector with which to compute the vector product
    /// \Returns The complex vector product of the two vectors
    ComplexVector3 Cross(ComplexVector3 vec);

    /// Real component
    TVector3 Real();

    /// Imaginary component
    TVector3 Imag();

    /// Return the complex conjugate of the vector
    ComplexVector3 Conj();

    /// Addition of vectors
    ComplexVector3 operator + (ComplexVector3 const &vec);

    /// Subtraction of vectors
    ComplexVector3 operator - (ComplexVector3 const &vec);

    /// Unary minus
    ComplexVector3 operator - () const;

    /// Multiplication with a real number
    ComplexVector3 operator * (double a);

    /// Multiplication with a complex number
    ComplexVector3 operator * (std::complex<double> a);
    
    /// Multiplication with a real number
    ComplexVector3 operator *= (double a);

    /// Multiplication with a complex number
    ComplexVector3 operator *= (std::complex<double> a);

  };
}

#endif
