#include "FFTWComplex.h"

#include <iostream>

std::ostream& operator<<(std::ostream& os, const rad::FFTWComplex& val) {
  return (os << "(" << val.re << "," << val.im << "i)");
}

double rad::getAbs(const rad::FFTWComplex& c) {
  return sqrt(c.re * c.re + c.im * c.im);
}