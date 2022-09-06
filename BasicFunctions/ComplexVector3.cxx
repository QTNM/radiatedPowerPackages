/// ComplexVector3.cxx

#include "BasicFunctions/ComplexVector3.h"

#include <complex>

#include "TVector3.h"

rad::ComplexVector3::ComplexVector3(TVector3 vec)
{
  mX = std::complex<double>{vec.X(), 0};
  mY = std::complex<double>{vec.Y(), 0};
  mZ = std::complex<double>{vec.Z(), 0};
}

std::complex<double> rad::ComplexVector3::Dot(ComplexVector3 vec)
{
  std::complex<double> sum{ mX*vec.X() + mY*vec.Y() + mZ*vec.Z() };
  return sum;
}

rad::ComplexVector3 rad::ComplexVector3::Cross(ComplexVector3 vec)
{
  ComplexVector3 product{ mY*vec.Z()-vec.Y()*mZ, mZ*vec.X()-vec.Z()*mX, mX*vec.Y()-vec.X()*mY };
  return product;
}

TVector3 rad::ComplexVector3::Real()
{
  return TVector3{ mX.real(), mY.real(), mZ.real() };
}

TVector3 rad::ComplexVector3::Imag()
{
  return TVector3{ mX.imag(), mY.imag(), mZ.imag() };
}

rad::ComplexVector3 rad::ComplexVector3::Conj()
{
  return ComplexVector3{ std::conj(mX), std::conj(mY), std::conj(mZ) };
}

rad::ComplexVector3 rad::ComplexVector3::operator + (ComplexVector3 const &vec)
{
  return ComplexVector3{ mX+vec.mX, mY+vec.mY, mZ+vec.mZ };
}

rad::ComplexVector3 rad::ComplexVector3::operator - (ComplexVector3 const &vec)
{
  return ComplexVector3{ mX-vec.mX, mY-vec.mY, mZ-vec.mZ };
}

rad::ComplexVector3 rad::ComplexVector3::operator - () const
{
  return ComplexVector3{ -mX, -mY, -mZ };
}

rad::ComplexVector3 rad::ComplexVector3::operator * (double a)
{
  return ComplexVector3{ a*mX, a*mY, a*mZ };
}

rad::ComplexVector3 rad::ComplexVector3::operator * (std::complex<double> a)
{
  return ComplexVector3{ a*mX, a*mY, a*mZ };
}

rad::ComplexVector3& rad::ComplexVector3::operator *= (double a)
{
  mX = mX * a;
  mY = mY * a;
  mZ = mZ * a;
  return *this;
}

rad::ComplexVector3& rad::ComplexVector3::operator *= (std::complex<double> a)
{
  mX = mX * a;
  mY = mY * a;
  mZ = mZ * a;
  return *this;
}
