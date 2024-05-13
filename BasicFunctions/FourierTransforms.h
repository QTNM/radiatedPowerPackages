/*
 * FourierTransforms.h
 *
 *  Created on: 13/05/2024
 */

#ifndef FOURIER_TRANSFORMS_H
#define FOURIER_TRANSFORMS_H

#include <fftw3.h>

#include "BasicFunctions/FFTWComplex.h"

namespace rad {
/// @brief Computes an FFT of an array of real numbers
/// @param length The length of the input array
/// @param theInput The input array of *length* real numbers
/// @return An array of *length/2 + 1* complex numbers
FFTWComplex *doFFT(int length, double *theInput);

/// @brief Computes an inverse FFT of an array of complex numbers
/// @param length The length of the output array
/// @param theInput The input array of complex number of *(length/2 + 1)*
/// @return An array of *length* real numbers
double *doInverseFFT(int length, const FFTWComplex *theInput);
}  // namespace rad

#endif  // FOURIER_TRANSFORMS_H