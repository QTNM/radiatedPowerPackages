/*
  ComplexVector3.h

  Header-only class supporting the use of complex three vector
  No external dependencies - pure C++ implementation
*/

#ifndef BASIC_CORE_COMPLEX_VECTOR_THREE_H
#define BASIC_CORE_COMPLEX_VECTOR_THREE_H

#include <array>
#include <complex>

namespace rad {

class ComplexVector3 {
 private:
  // The components
  std::complex<double> mX, mY, mZ;

 public:
  /// Parametrised constructor taking three complex numbers
  inline ComplexVector3(std::complex<double> x, std::complex<double> y,
                        std::complex<double> z)
      : mX{x}, mY{y}, mZ{z} {}

  /// Parametrised constructor taking six real numbers (real and imaginary
  /// parts)
  inline ComplexVector3(double xReal, double yReal, double zReal,
                        double xImag = 0.0, double yImag = 0.0,
                        double zImag = 0.0)
      : mX{xReal, xImag}, mY{yReal, yImag}, mZ{zReal, zImag} {}

  /// Functions to return individual components
  inline std::complex<double> X() const { return mX; }
  inline std::complex<double> Y() const { return mY; }
  inline std::complex<double> Z() const { return mZ; }

  /// Functions to set components after creation
  inline void SetX(std::complex<double> x) { mX = x; }
  inline void SetY(std::complex<double> y) { mY = y; }
  inline void SetZ(std::complex<double> z) { mZ = z; }

  /// Get real components as array
  inline std::array<double, 3> RealComponents() const {
    return {mX.real(), mY.real(), mZ.real()};
  }

  /// Get imaginary components as array
  inline std::array<double, 3> ImagComponents() const {
    return {mX.imag(), mY.imag(), mZ.imag()};
  }

  /// Scalar product
  /// \param vec The vector with which to compute the scalar product
  /// \Returns The complex scalar product of the two vectors
  inline std::complex<double> Dot(const ComplexVector3& vec) const {
    return mX * vec.mX + mY * vec.mY + mZ * vec.mZ;
  }

  /// Vector product
  /// \param vec The vector with which to compute the vector product
  /// \Returns The complex vector product of the two vectors
  inline ComplexVector3 Cross(const ComplexVector3& vec) const {
    return ComplexVector3{mY * vec.mZ - vec.mY * mZ, mZ * vec.mX - vec.mZ * mX,
                          mX * vec.mY - vec.mX * mY};
  }

  /// Return the complex conjugate of the vector
  inline ComplexVector3 Conj() const {
    return ComplexVector3{std::conj(mX), std::conj(mY), std::conj(mZ)};
  }

  /// Addition of vectors
  inline ComplexVector3 operator+(const ComplexVector3& vec) const {
    return ComplexVector3{mX + vec.mX, mY + vec.mY, mZ + vec.mZ};
  }

  /// Subtraction of vectors
  inline ComplexVector3 operator-(const ComplexVector3& vec) const {
    return ComplexVector3{mX - vec.mX, mY - vec.mY, mZ - vec.mZ};
  }

  /// Unary minus
  inline ComplexVector3 operator-() const {
    return ComplexVector3{-mX, -mY, -mZ};
  }

  /// Multiplication with a real number
  inline ComplexVector3 operator*(double a) const {
    return ComplexVector3{mX * a, mY * a, mZ * a};
  }

  /// Multiplication with a complex number
  inline ComplexVector3 operator*(std::complex<double> a) const {
    return ComplexVector3{mX * a, mY * a, mZ * a};
  }

  /// Multiplication with a real number
  inline ComplexVector3& operator*=(double a) {
    mX *= a;
    mY *= a;
    mZ *= a;
    return *this;
  }

  /// Multiplication with a complex number
  inline ComplexVector3& operator*=(std::complex<double> a) {
    mX *= a;
    mY *= a;
    mZ *= a;
    return *this;
  }

  /// Compound assignment addition
  inline ComplexVector3& operator+=(const ComplexVector3& vec) {
    mX += vec.mX;
    mY += vec.mY;
    mZ += vec.mZ;
    return *this;
  }

  /// Compound assignment subtraction
  inline ComplexVector3& operator-=(const ComplexVector3& vec) {
    mX -= vec.mX;
    mY -= vec.mY;
    mZ -= vec.mZ;
    return *this;
  }
};

/// Multiplication by a double (left side)
inline ComplexVector3 operator*(double a, const ComplexVector3& vec) {
  return vec * a;
}

/// Multiplication by a complex double (left side)
inline ComplexVector3 operator*(const std::complex<double>& a,
                                const ComplexVector3& vec) {
  return vec * a;
}

}  // namespace rad

#endif